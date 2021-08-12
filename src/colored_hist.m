function varargout = colored_hist(data,bins)
%COLORED_HIST plots a histogram and color codes the bins according to their
%value.

% Histogram
[n,c] = hist(data,bins);
freq = n/sum(n);
% if max(freq) > ymax
%     ymax = max(freq);
% end
% Make faces and vertices for plotting the histograms
faces = NaN(length(n),4);
vertices = zeros(length(n)*4,2);
w = c(2)-c(1);
for i = 1 : length(c)
    vertices((1:4) + 4*(i-1),:) = [[1 -1 -1 1]'*w/2+c(i) [0 0 1 1]'*freq(i)];
    faces(i,:) = [1 2 3 4] + (i-1)*4;
end
h =patch('Vertices',vertices,...
    'Faces',faces,...
    'FaceColor','flat',...
    'FaceVertexCData',c',...
    'EdgeColor','k',...
    'EdgeAlpha',1,...
    'LineWidth',1);

if nargout>0
    varargout{1} = h;
end

end

