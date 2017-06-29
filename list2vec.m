function [countvec,timefirevec] = list2vec(list,timevec,varargin)

% DR - added a varargin marker to input an array the length of list
% whose number at each row marker(x) gets input into the output array for all the
% vec indices between list(x, [start end]) times
if (~isempty(varargin))
    assign(varargin{:});
end

% vector = list2vec(list,timevec)
% converts list spiketimes in timevec
% into vector of 0-1 for corresponding time bins in list

% initialize
countvec = zeros(length(list),1);
timefirevec = zeros(length(list),1);

for p = 1:size(timevec)
    startind = find(list > timevec(p), 1, 'first') - 1;
    countvec(startind) = countvec(startind) + 1;
    timefirevec(startind) = timevec(p);
end

end
