function JOBLIST = PIVJobList()

% Data: Input images
Data.Inputs.Images.Directory = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/pivchallenge/2014/A/images/proc/ghost';
Data.Inputs.Images.BaseName = 'A_deghost_';
Data.Inputs.Images.Digits = 5;
Data.Inputs.Images.Extension = '.tif';
Data.Inputs.Images.Trailers = {'_a', '_b'};
Data.Inputs.Images.Frames.Start = 1;
Data.Inputs.Images.Frames.End = 300;
Data.Inputs.Images.Frames.Step = 1;

% Data: output vectors
Data.Outputs.Vectors.Directory = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/pivchallenge/2014/A/vect';
Data.Outputs.Vectors.BaseName = 'A_deghost_';
Data.Outputs.Vectors.Digits = 5;
Data.Outputs.Vectors.Extension = '.mat';

% Interrogation region dimensions
Processing(1).Region.Height = 64;
Processing(1).Region.Width = 128;

% Spatial window
Processing(1).Window.Fraction = {[0.25, 0.5], [0.25, 1]};

% Grid parameters
Processing(1).Grid.Spacing.Y = 64;
Processing(1).Grid.Spacing.X = 64;
Processing(1).Grid.Shift.Y = -16;
Processing(1).Grid.Shift.X = 0;
Processing(1).Grid.Mask.Directory = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/pivchallenge/2014/A/images/masks';
Processing(1).Grid.Mask.Name = 'imgAmask3.tif';

% Correlation parameters
Processing(1).Correlation.EstimatedParticleDiameter = 6;
Processing(1).Correlation.Method = 'apc';
Processing(1).Correlation.Step = 0;
Processing(1).Correlation.Ensemble.NumberOfPairs = 10;

% Parameters specific to APC
Processing(1).Correlation.APC.EnsembleLength = 10;
Processing(1).Correlation.APC.Shuffle.Range = [0, 0];
Processing(1).Correlation.APC.Shuffle.Step = [0, 0];

% Parameters specific to RPC
Processing(1).Correlation.RPC.FilterDiameter = 6;

% Parameters for vector validation
Processing(1).Validation.DoValidation = 1;
Processing(1).Validation.ValidationMethod = 'uod';

% Parameters for smoothing
Processing(1).Smoothing.DoSmoothing = 1;
Processing(1).Smoothing.KernelDiameter = 7;
Processing(1).Smoothing.KernelStdDev = 1;

% Parameters for iterative method
Processing(1).Iterative.Method = 'deform';
Processing(1).Iterative.Deform.Interpolation = 'interp2'; 
Processing(1).Iterative.Deform.ConvergenceCriterion = 0.1;
Processing(1).Iterative.Deform.MaxIterations = 1;

% Copy the first processing pass to the second one.
Processing(2) = Processing(1);

% Add the fields to the jobfile structure.
JobFile.Data = Data;
JobFile.Processing = Processing;

JOBLIST(1) = JobFile;

end



