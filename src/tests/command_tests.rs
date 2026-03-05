#[cfg(test)]
mod tests {
    use std::io::Write;
    use tempfile::NamedTempFile;
    use rstest::{fixture, rstest};

    use crate::{
        errors::PisError,
        readers::{
            input_file::commands::Command,
            simulation_context::{SimulationContext, VelocityDistribution},
        },
    };

    // ─── fixtures ───────────────────────────────────────────────────────────────

    #[fixture]
    fn ctx() -> SimulationContext {
        SimulationContext::default()
    }

    /// A context that already has pair_style initialised — used by pair_coeff tests.
    #[fixture]
    fn ctx_with_pair_style(mut ctx: SimulationContext) -> SimulationContext {
        run_cmd_ok(&Command::PairStyle, "lj/cut 10.0", &mut ctx);
        ctx
    }

    // ─── helpers ────────────────────────────────────────────────────────────────

    fn run_cmd_ok(cmd: &Command, raw_args: &str, ctx: &mut SimulationContext) {
        let args: Vec<&str> = raw_args.split_whitespace().collect();
        cmd.run(&args, 0, ctx).expect("expected command to succeed");
    }

    fn run_cmd_err(cmd: &Command, raw_args: &str, ctx: &mut SimulationContext) -> PisError {
        let args: Vec<&str> = raw_args.split_whitespace().collect();
        cmd.run(&args, 0, ctx)
            .expect_err("expected command to return an error, but it succeeded")
    }

    fn write_temp_data_file(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("could not create temp file");
        write!(f, "{}", content).unwrap();
        f
    }

    // ─── Command::from_str ──────────────────────────────────────────────────────

    #[rstest]
    #[case("timestep")]
    #[case("run")]
    #[case("velocity")]
    #[case("read_data")]
    #[case("fix")]
    #[case("pair_style")]
    #[case("pair_coeff")]
    #[case("dump")]
    fn from_str_known_commands_return_some(#[case] input: &str) {
        assert!(Command::from_str(input).is_some(), "expected Some for '{input}'");
    }

    #[rstest]
    #[case("unknown_command")]
    #[case("")]
    #[case("TIMESTEP")]
    #[case("TimeStep")]
    fn from_str_unknown_command_returns_none(#[case] input: &str) {
        assert!(Command::from_str(input).is_none(), "expected None for '{input}'");
    }

    // ─── timestep ───────────────────────────────────────────────────────────────

    #[rstest]
    fn timestep_sets_value_on_ctx(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::TimeStep, "0.001", &mut ctx);
        assert!((ctx.timestep - 0.001).abs() < f64::EPSILON);
    }

    #[rstest]
    fn timestep_missing_arg_returns_missing_argument(mut ctx: SimulationContext) {
        let err = run_cmd_err(&Command::TimeStep, "", &mut ctx);
        assert!(
            matches!(&err, PisError::MissingArgument { line: 0 }),
            "expected MissingArgument {{ line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    #[case("notanumber")]
    #[case("one")]
    #[case("0.1.2")]
    fn timestep_invalid_float_returns_float_parse_error(mut ctx: SimulationContext, #[case] input: &str) {
        let err = run_cmd_err(&Command::TimeStep, input, &mut ctx);
        assert!(
            matches!(&err, PisError::FloatParseError { string, .. } if string == input),
            "expected FloatParseError for \"{input}\", got: {err:?}"
        );
    }

    // ─── run (RunSteps) ─────────────────────────────────────────────────────────

    #[rstest]
    fn runsteps_sets_value_on_ctx(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::RunSteps, "5000", &mut ctx);
        assert_eq!(ctx.steps, 5000);
    }

    #[rstest]
    fn runsteps_missing_arg_returns_missing_argument(mut ctx: SimulationContext) {
        let err = run_cmd_err(&Command::RunSteps, "", &mut ctx);
        assert!(
            matches!(&err, PisError::MissingArgument { line: 0 }),
            "expected MissingArgument {{ line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    fn runsteps_negative_value_returns_negative_value_error(mut ctx: SimulationContext) {
        let err = run_cmd_err(&Command::RunSteps, "-1", &mut ctx);
        assert!(
            matches!(&err, PisError::NegativeValue { value: -1, line: 0 }),
            "expected NegativeValue {{ value: -1, line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    #[case("1.5")]
    #[case("ten")]
    #[case("1e3")]
    fn runsteps_non_integer_returns_int_parse_error(mut ctx: SimulationContext, #[case] input: &str) {
        let err = run_cmd_err(&Command::RunSteps, input, &mut ctx);
        assert!(
            matches!(&err, PisError::IntParseError { string, .. } if string == input),
            "expected IntParseError for \"{input}\", got: {err:?}"
        );
    }

    // ─── velocity ───────────────────────────────────────────────────────────────

    #[rstest]
    fn velocity_create_minimal_sets_ctx(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::Velocity, "all create 300.0 12345", &mut ctx);
        let sv = ctx.starting_velocity.expect("starting_velocity should be Some");
        assert_eq!(sv.group, "all");
        assert!((sv.start_temperature.unwrap() - 300.0).abs() < f64::EPSILON);
        assert_eq!(sv.seed.unwrap(), 12345);
    }

    #[rstest]
    fn velocity_create_minimal_sets_without_seed_ctx(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::Velocity, "all create 300.0", &mut ctx);
        let sv = ctx.starting_velocity.expect("starting_velocity should be Some");
        assert_eq!(sv.group, "all");
        assert!((sv.start_temperature.unwrap() - 300.0).abs() < f64::EPSILON);
        assert_eq!(sv.seed.unwrap(), 0);
    }

    #[rstest]
    #[case("all create 300.0 42 dist gaussian", VelocityDistribution::Gaussian)]
    #[case("all create 300.0 42 dist uniform",  VelocityDistribution::Uniform)]
    fn velocity_create_dist_keyword_sets_distribution(
        mut ctx: SimulationContext,
        #[case] args: &str,
        #[case] expected_dist: VelocityDistribution,
    ) {
        run_cmd_ok(&Command::Velocity, args, &mut ctx);
        let sv = ctx.starting_velocity.unwrap();
        assert_eq!(sv.dist, Some(expected_dist));
    }

    #[rstest]
    fn velocity_create_dist_keyword_sets_distribution_without_seed(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::Velocity, "all create 300.0 dist gaussian", &mut ctx);
        let sv = ctx.starting_velocity.unwrap();
        assert_eq!(sv.dist, Some(VelocityDistribution::Gaussian));
    }

    #[rstest]
    fn velocity_missing_group_returns_missing_argument(mut ctx: SimulationContext) {
        let err = run_cmd_err(&Command::Velocity, "", &mut ctx);
        assert!(
            matches!(&err, PisError::MissingArgument { line: 0 }),
            "expected MissingArgument {{ line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    #[case("all set 300.0 42",   "set")]
    #[case("all scale 300.0 42", "scale")]
    fn velocity_unknown_style_returns_invalid_argument(
        mut ctx: SimulationContext,
        #[case] args: &str,
        #[case] bad_token: &str,
    ) {
        let err = run_cmd_err(&Command::Velocity, args, &mut ctx);
        assert!(
            matches!(&err, PisError::InvalidArgument { string, line: 0 } if string == bad_token),
            "expected InvalidArgument {{ string: \"{bad_token}\", line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    #[case("all create 300.0 42 badkeyword value", "badkeyword")]
    #[case("all create 300.0 42 dist lorentzian",  "lorentzian")]
    fn velocity_invalid_keyword_or_dist_returns_invalid_argument(
        mut ctx: SimulationContext,
        #[case] args: &str,
        #[case] bad_token: &str,
    ) {
        let err = run_cmd_err(&Command::Velocity, args, &mut ctx);
        assert!(
            matches!(&err, PisError::InvalidArgument { string, line: 0 } if string == bad_token),
            "expected InvalidArgument {{ string: \"{bad_token}\", line: 0 }}, got: {err:?}"
        );
    }

    // ─── pair_style ─────────────────────────────────────────────────────────────

    #[rstest]
    fn pair_style_stores_args_in_ctx(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::PairStyle, "lj/cut 10.0", &mut ctx);
        let pa = ctx.potential_args.expect("potential_args should be Some");
        assert_eq!(pa.pair_style_args, vec!["lj/cut", "10.0"]);
    }

    #[rstest]
    fn pair_style_records_line_number(mut ctx: SimulationContext) {
        let args: Vec<&str> = vec!["lj/cut", "10.0"];
        Command::PairStyle.run(&args, 7, &mut ctx).unwrap();
        assert_eq!(ctx.potential_args.unwrap().pair_style_line, 7);
    }

    // ─── pair_coeff ─────────────────────────────────────────────────────────────

    #[rstest]
    fn pair_coeff_without_prior_pair_style_returns_potential_not_initialized(mut ctx: SimulationContext) {
        let err = run_cmd_err(&Command::PairCoeff, "1 1 0.238 3.405 8.5", &mut ctx);
        assert!(
            matches!(&err, PisError::PotentialNotInitialized),
            "expected PotentialNotInitialized, got: {err:?}"
        );
    }

    #[rstest]
    fn pair_coeff_appends_to_existing_potential_args(mut ctx_with_pair_style: SimulationContext) {
        run_cmd_ok(&Command::PairCoeff, "1 1 0.238 3.405 8.5", &mut ctx_with_pair_style);
        run_cmd_ok(&Command::PairCoeff, "1 2 0.150 3.100 7.5", &mut ctx_with_pair_style);
        let pa = ctx_with_pair_style.potential_args.unwrap();
        assert_eq!(pa.pair_coeff_args.len(), 2);
        assert_eq!(pa.pair_coeff_args[0], vec!["1", "1", "0.238", "3.405", "8.5"]);
        assert_eq!(pa.pair_coeff_args[1], vec!["1", "2", "0.150", "3.100", "7.5"]);
    }

    // ─── dump ───────────────────────────────────────────────────────────────────

    #[rstest]
    fn dump_sets_all_fields_on_ctx(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::Dump, "my_dump all atom 100 output.lammpstrj", &mut ctx);
        assert_eq!(ctx.dump_args.name, "my_dump");
        assert_eq!(ctx.dump_args.group, "all");
        assert_eq!(ctx.dump_args.style, "atom");
        assert_eq!(ctx.dump_args.dump_step, 100);
        assert_eq!(ctx.dump_args.file_name, "output.lammpstrj");
    }

    #[rstest]
    #[case("my_dump all atom 100")]  // missing file_name
    #[case("my_dump all atom")]      // missing dump_step and file_name
    #[case("my_dump all")]           // missing style, dump_step, and file_name
    fn dump_missing_args_returns_missing_argument(mut ctx: SimulationContext, #[case] args: &str) {
        let err = run_cmd_err(&Command::Dump, args, &mut ctx);
        assert!(
            matches!(&err, PisError::MissingArgument { line: 0 }),
            "expected MissingArgument {{ line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    #[case("my_dump all atom every output.lammpstrj", "every")]
    #[case("my_dump all atom 1.5 output.lammpstrj",   "1.5")]
    fn dump_non_integer_step_returns_int_parse_error(
        mut ctx: SimulationContext,
        #[case] args: &str,
        #[case] bad_token: &str,
    ) {
        let err = run_cmd_err(&Command::Dump, args, &mut ctx);
        assert!(
            matches!(&err, PisError::IntParseError { string, .. } if string == bad_token),
            "expected IntParseError for \"{bad_token}\", got: {err:?}"
        );
    }

    // ─── fix ────────────────────────────────────────────────────────────────────

    #[rstest]
    fn fix_nvt_temp_keyword_sets_nh_chain_args(mut ctx: SimulationContext) {
        run_cmd_ok(&Command::Fix, "myfix all nvt temp 300.0 300.0 0.1", &mut ctx);
        let nh = ctx.nh_chain_args.expect("nh_chain_args should be set");
        assert_eq!(nh.name, "myfix");
        assert_eq!(nh.group, "all");
        assert!((nh.start_temperature - 300.0).abs() < f64::EPSILON);
        assert!((nh.tau - 0.1).abs() < f64::EPSILON);
    }

    #[rstest]
    fn fix_npt_sets_both_barostat_and_thermostat(mut ctx: SimulationContext) {
        run_cmd_ok(
            &Command::Fix,
            "myfix all npt temp 300.0 300.0 0.1 iso 1.0 1.0 1.0",
            &mut ctx,
        );
        assert!(ctx.nh_chain_args.is_some(), "thermostat args should be set");
        assert!(ctx.mtk_barostat_args.is_some(), "barostat args should be set");
    }

    #[rstest]
    fn fix_npt_iso_sets_barostat_pressure_and_tau(mut ctx: SimulationContext) {
        run_cmd_ok(
            &Command::Fix,
            "myfix all npt temp 300.0 300.0 0.1 iso 2.0 2.0 0.5",
            &mut ctx,
        );
        let baro = ctx.mtk_barostat_args.unwrap();
        assert!((baro.start_pressure[(0, 0)] - 2.0).abs() < f64::EPSILON);
        assert!((baro.tau - 0.5).abs() < f64::EPSILON);
    }

    #[rstest]
    #[case("myfix all")]  // style missing
    #[case("myfix")]      // group and style missing
    fn fix_missing_required_args_returns_missing_argument(
        mut ctx: SimulationContext,
        #[case] args: &str,
    ) {
        let err = run_cmd_err(&Command::Fix, args, &mut ctx);
        assert!(
            matches!(&err, PisError::MissingArgument { line: 0 }),
            "expected MissingArgument {{ line: 0 }}, got: {err:?}"
        );
    }

    #[rstest]
    #[case("myfix all nvt badkeyword 300.0 300.0 0.1", "badkeyword")]
    #[case("myfix all npt unknown 1.0 1.0 1.0",        "unknown")]
    fn fix_unknown_keyword_returns_invalid_argument(
        mut ctx: SimulationContext,
        #[case] args: &str,
        #[case] bad_token: &str,
    ) {
        let err = run_cmd_err(&Command::Fix, args, &mut ctx);
        assert!(
            matches!(&err, PisError::InvalidArgument { string, line: 0 } if string == bad_token),
            "expected InvalidArgument {{ string: \"{bad_token}\", line: 0 }}, got: {err:?}"
        );
    }

    // ─── read_data ──────────────────────────────────────────────────────────────

    #[rstest]
    fn read_data_populates_atoms_and_box(mut ctx: SimulationContext) {
        let content = "\
# Simple test data file

3 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 1.0

Atoms

1 1 1.0 2.0 3.0
2 1 4.0 5.0 6.0
3 1 7.0 8.0 9.0
";
        let f = write_temp_data_file(content);
        let path = f.path().to_str().unwrap();
        Command::ReadData.run(&[path], 0, &mut ctx).unwrap();
        let atoms = ctx.atoms.expect("atoms should be Some");
        assert_eq!(atoms.n_atoms, 3);
        assert!((atoms.positions[(0, 0)] - 1.0).abs() < f64::EPSILON);
        assert!((atoms.positions[(1, 1)] - 5.0).abs() < f64::EPSILON);
    }

    #[rstest]
    fn read_data_nonexistent_file_returns_input_file_error(mut ctx: SimulationContext) {
        let bad_path = "/nonexistent/path/to/file.data";
        let err = Command::ReadData
            .run(&[bad_path], 0, &mut ctx)
            .expect_err("expected an error for a missing file");
        assert!(
            matches!(&err, PisError::InputFileError { path, .. } if path == bad_path),
            "expected InputFileError for path \"{bad_path}\", got: {err:?}"
        );
    }

    #[rstest]
    #[case(0, 2)]  // atom ID 0 is invalid (IDs are 1-indexed)
    #[case(5, 2)]  // atom ID 5 is out of range for a 2-atom system
    fn read_data_invalid_atom_id_returns_atom_count_mismatch(
        mut ctx: SimulationContext,
        #[case] bad_id: usize,
        #[case] n_atoms: usize,
    ) {
        let content = format!(
            "\
{n_atoms} atoms
1 atom types

0.0 5.0 xlo xhi
0.0 5.0 ylo yhi
0.0 5.0 zlo zhi

Masses

1 1.0

Atoms

1 1 0.0 0.0 0.0
{bad_id} 1 1.0 1.0 1.0
"
        );
        let f = write_temp_data_file(&content);
        let path = f.path().to_str().unwrap();
        let err = Command::ReadData
            .run(&[path], 0, &mut ctx)
            .expect_err("expected an error for an invalid atom ID");
        assert!(
            matches!(&err, PisError::AtomCountMismatch { expected, found } if *expected == n_atoms && *found == bad_id),
            "expected AtomCountMismatch {{ expected: {n_atoms}, found: {bad_id} }}, got: {err:?}"
        );
    }

    #[rstest]
    fn read_data_with_velocities_sets_start_velocity_flag_false(mut ctx: SimulationContext) {
        let content = "\
2 atoms
1 atom types

0.0 5.0 xlo xhi
0.0 5.0 ylo yhi
0.0 5.0 zlo zhi

Masses

1 1.0

Atoms

1 1 0.0 0.0 0.0
2 1 1.0 1.0 1.0

Velocities

1 0.1 0.2 0.3
2 0.4 0.5 0.6
";
        let f = write_temp_data_file(content);
        let path = f.path().to_str().unwrap();
        Command::ReadData.run(&[path], 0, &mut ctx).unwrap();
        let sv = ctx.starting_velocity.expect("starting_velocity should be Some");
        assert!(
            !sv.start_velocity,
            "start_velocity should be false when a Velocities section is present"
        );
    }

    #[rstest]
    fn read_data_with_pair_coeffs_populates_potential_manager(mut ctx: SimulationContext) {
        let content = "\
2 atoms
1 atom types

0.0 5.0 xlo xhi
0.0 5.0 ylo yhi
0.0 5.0 zlo zhi

Masses

1 1.0

PairCoeffs

1 0.238 3.405 8.5

Atoms

1 1 0.0 0.0 0.0
2 1 1.0 1.0 1.0
";
        let f = write_temp_data_file(content);
        let path = f.path().to_str().unwrap();
        Command::ReadData.run(&[path], 0, &mut ctx).unwrap();
        assert!(ctx.mgr.is_some(), "potential manager should be populated from PairCoeffs");
    }
}
