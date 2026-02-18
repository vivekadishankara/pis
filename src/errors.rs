use thiserror::Error;

#[allow(dead_code)]
#[derive(Error, Debug)]
pub enum PisError {
    // File I/O Errors
    #[error("Failed to open input file '{path}': {source})")]
    InputFileError {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("Failed to read line {line} in file '{path}': {source}")]
    DataFileError {
        path: String,
        line: usize,
        #[source]
        source: std::io::Error,
    },

    // Parsing Errors
    #[error("Invalid command {command} found line: {line}")]
    UnkownCommand {
        command: String,
        line: usize,
    },

    #[error("Missing argument on line {line}")]
    MissingArgument {line: usize },

    #[error("Error parsing floating number from string {string}: {source}")]
    FloatParseError {
        string: String,
        #[source]
        source: std::num::ParseFloatError,
    },

    #[error("Error parsing integer number from string {string}: {source}")]
    IntParseError {
        string: String,
        #[source]
        source: std::num::ParseIntError,
    },

    #[error("Negative value {value} not allowed on line: {line}")]
    NegativeValue {
        value: i32,
        line: usize,
    },

    #[error("Invalid argument: {string} at line: {line}")]
    InvalidArgument { string: String, line: usize },

    // Configuration errors
    #[error("No atoms defined in input file")]
    NoAtomsDefined,

    #[error("Atom count mismatch: expected {expected}, found {found}")]
    AtomCountMismatch { expected: usize, found: usize },

    #[error("Potential manager not initialized - missing pair_style or pair_coeff commands")]
    PotentialNotInitialized,

    #[error("Unknown pair style: '{style}'")]
    UnknownPairStyle { style: String },

    #[error("Missing potential for atom pair ({i}, {j})")]
    MissingPotential { i: usize, j: usize },

    // Physics errors
    #[error("Atom {id} has NaN enrgy at step {step}")]
    NaNEnergy { id: usize, step: usize },

    // Array bounds errors
    #[error("Atom type {type_id} out of range")]
    InvalidAtomType { type_id: usize },

    #[error("Atom index {index} out of range (total atoms: {n_atoms})")]
    InvalidAtomIndex { index: usize, n_atoms: usize },

    // #[error(transparent)]
    // Other(#[from] anyhow::Error),
}

pub type Result<T> = std::result::Result<T, PisError>;
