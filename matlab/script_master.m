addpath('half-edge-mesh');

warning('off','MATLAB:MKDIR:DirectoryExists');

mkdir('pavilion');
test_pavilion;

mkdir('paper-spiral');
test_paper_spiral_gravity;

test_revolution_examples;

test_horse_0;
test_horse_1;
test_horse_2;
test_horse_3;
plot_horse_figure;

test_flower_pot_1;
test_flower_pot_2;
plot_flower_pot;
plot_gravity_refinement;

plot_stability_example;
plot_stability_opt;

mkdir('ellipse-figures');
plot_ellipse_example;

test_teaser_1;
test_teaser_2;
plot_teaser;

mkdir('paper-stability-S');
test_stability_S;
mkdir('robustness');
test_robustness_1;
test_robustness_2;