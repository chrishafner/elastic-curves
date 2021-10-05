addpath('half-edge-mesh');

warning('off','MATLAB:MKDIR:DirectoryExists');

% Demo 1: Pavilion. Reproduces Fig. 17, and shows a rendering of the
% model that serves as the basis for the physical model shown in Fig. 1
% (left). An SVG file with the cutting paths for the cardboard model is
% saved to the 'pavilion/strips.svg'. In addition, the script generates
% 'pavilion/pavilion_base.txt', containing FeatureScript code to generate
% proxy geometry for the slots in the support structure in OnShape.
mkdir('pavilion');
test_pavilion;

% Demo 2: Spiral. Reproduces Fig. 7 and stores the cutting paths for the
% models with and without taking gravity into account to the subfolder
% 'paper-spiral'.
mkdir('paper-spiral');
test_paper_spiral_gravity;

% Demo 3: Surfaces of Revolution. Reproduces Fig. 22.
test_revolution_examples;

% Demo 4: Horse. Reproduces Fig. 18, and a rendering of the model that
% serves as the basis for the physical model shown in Fig. 24 (top right).
% Also stores the cutting paths and FeatureScript for proxy geometry to the
% 'horse' subfolder.
test_horse_0;
test_horse_1;
test_horse_2;
test_horse_3;
plot_horse_figure;

% Demo 5: Flower Pot. Reproduces Figs. 19 and 8, and stores cutting paths
% to the 'flower-pot' subfolder.
test_flower_pot_1;
test_flower_pot_2;
plot_flower_pot;
plot_gravity_refinement;

% Demo 6: Stability. Reproduces Figs. 9 and 10.
plot_stability_example;
plot_stability_opt;

% Demo 7: Stiffness Family. Reprduces Fig. 6.
mkdir('ellipse-figures');
plot_ellipse_example;

% Demo 8: Shell. Reproduces Fig. 20, and shows cutting paths for the
% cardboard model.  The result is stored to the 'teaser' subfolder.
test_teaser_1;
test_teaser_2;
plot_teaser;

% Demo 9: Robustness Tests. Reprocuces Fig. 15.
mkdir('paper-stability-S');
test_stability_S;
mkdir('robustness');
test_robustness_1;
test_robustness_2;