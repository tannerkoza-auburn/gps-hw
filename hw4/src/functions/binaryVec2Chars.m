function [charVec] = binaryVec2Chars(binaryVec)
%BINARYVEC2CHARS Accepts a binary Vector
% The function will rashape it into a 8xN vector
% and output the string message and return the message vector

outputVec = reshape(binaryVec,[8 11]);
charVec = [];
for i = 1:length(outputVec)
    str = strjoin(string(outputVec(:,i)));
    decVal = bin2dec(str);
    charVal = char(decVal);
    charVec(i) = charVal;
end
fprintf('The Message Says: %s\n',char(charVec));
end

