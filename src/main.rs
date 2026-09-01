mod mode;
mod node_math;
mod graph;
mod channel;
mod random_ldpc;
mod css;
mod decoder;

mod sim_classic;
mod sim_jd;
mod sim_css;

fn main() {
    //sim_classic::run();
    //sim_jd::run();
    sim_css::run();
}



