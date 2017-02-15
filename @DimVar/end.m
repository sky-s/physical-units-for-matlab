function endIndex = end(vIn,indexToFindEnd,numberOfIndexes)
sizeVector = size(vIn.value);

if(numberOfIndexes == 1)
    endIndex = max(sizeVector);
else
    endIndex = sizeVector(indexToFindEnd);
end
