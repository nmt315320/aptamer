function dHRC= dHRC(DNAx,DNAy )
%����DNAx��DNAy���油���еĺ�������,�Ա����Գƺ���
xLen=size(DNAx,2);
yLen=size(DNAy,2);
if xLen==yLen
    dHRC=size(find(abs(DNAx-(5-flip(DNAy)))>0),2);
else
    error('����DNA���Ȳ�ͬ');
end
end

