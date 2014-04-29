% the ID of the B-scan (for Spectralis scans either 0 for 2-D or {0,...,61} for volumes)
collector.options.labelID = 0;
collector.options.loadRoutineData = 'spectralis';
collector.options.folder_data = '/export/home/oct/OCT-labeled-testdata/mat/';
% if restore == 1 and we have already labeled data with the label tool
collector.options.saveDir = '/export/home/oct/OCT-labeled-testdata/GT-Mehrschicht/stefan_backup/'; 
% if restore == 1, we have not labeled anything yet but there are labels provided with the dataset
collector.options.loadRoutineLabels = 'spectralisLabels'; 
collector.options.folder_labels = '/export/home/oct/OCT-labeled-testdata/GT-Mehrschicht/stefan_backup/';
% we want no preprocessing (we could add here a Gaussian filter for example or any user-defined function that eases labeling)
collector.options.preprocessing.scanLevel = {};
% using existing labels?
restore = 1;
% we label the B-Scan with ID stored inside 'filename.mat';
labelTool('filename',collector,restore);
