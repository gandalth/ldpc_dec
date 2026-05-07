use std::fmt;

pub enum OpMode {
    Classic,
    Quantum,
}

impl fmt::Display for OpMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OpMode::Classic => write!(f, "Classic Linear Code"),
            OpMode::Quantum => write!(f, "Quantum Error Correction (QEC) Code"),
        }
    }
}
