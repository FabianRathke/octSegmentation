function f1 = displayAppearance(B0,appearance,columnsPred)
	f1 = figure;
    h1 = imagesc(B0); colormap('gray'); hold on;
    app = zeros(size(B0));
    app(:,columnsPred) = appearance;
    % normalize the maximum of each column to one
%    app(:,columnsPred) = app(:,columnsPred)./repmat(max(app(:,columnsPred)),size(B0,1),1);
    % interpolate columns
%    columnsToInterpolate = setdiff(columnsPred(1):columnsPred(end),columnsPred);
%
%    [X Y] = meshgrid(columnsPred,1:size(B0,1));
%    [Xq Yq] = meshgrid(columnsToInterpolate,1:size(B0,1));
%    app(:,columnsToInterpolate) = interp2(X,Y,app(:,columnsPred),Xq,Yq);

    h2 = imagesc(app);
    set(h2,'AlphaData',app);
    colormap(f1,[gray(64); jet(64)])
    set(h1,'CData',get(h1,'CData')/2)
    set(h2,'CData',(get(h2,'CData')+1)/2)
end

