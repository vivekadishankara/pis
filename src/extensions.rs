use crate::errors::{PisError, Result};

pub trait ArgsExt {
    fn get_required(&self, index: usize, line: usize) -> Result<&str>;
    fn parse_int_at(&self, index: usize, line: usize) -> Result<i32>;
    fn parse_float_at(&self, index: usize, line: usize) -> Result<f64>;
}

impl ArgsExt for [&str] {
    fn get_required(&self, index: usize, line: usize) -> Result<&str> {
        self.get(index)
            .copied()
            .ok_or(PisError::MissingArgument { line })
    }
    
    fn parse_int_at(&self, index: usize, line: usize) -> Result<i32> {
        let arg = self.get_required(index, line)?;
        arg.parse()
            .map_err(|e| PisError::IntParseError {
                string: arg.to_string(),
                source: e,
            })
    }

    fn parse_float_at(&self, index: usize, line: usize) -> Result<f64> {
        let arg = self.get_required(index, line)?;
        arg.parse()
            .map_err(|e| PisError::FloatParseError { 
                string: arg.to_string(),
                source: e 
            })
    }
}

pub trait Int32ToUsize {
    fn convert_to_usize(&self, line: usize) -> Result<usize>;
}

impl Int32ToUsize for i32 {
    fn convert_to_usize(&self, line: usize) -> Result<usize> {
        (*self).try_into()
            .map_err(|_| PisError::NegativeValue{
                value: (*self),
                line: line
            })
    }
}