function B0 = addNoise(B0,options)

fprintf('Add speckle noise with variance %.2f to the image\n',options{2});
% A = A + n*A where n is a uniform distribution with mean 0 and variance options{2}
B0 = max(0,min(B0 + (rand(size(B0))-0.5)*sqrt(12*options{2}).*B0,1)); 
