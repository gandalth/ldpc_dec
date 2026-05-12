pub trait Mode {
    fn name() -> &'static str;
}

pub struct Classic;
pub struct Quantum;

impl Mode for Classic {
    fn name() -> &'static str {
        "Classic"
    }
}

impl Mode for Quantum {
    fn name() -> &'static str {
        "Quantum"
    }
}
