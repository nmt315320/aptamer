function distHammingReverse = dHR( DNAx,DNAy )
%����DNAx��DNAy�������еĺ�������,�Ա����Գƺ���
xLen=size(DNAx,2);
yLen=size(DNAy,2);
if xLen==yLen
    distHammingReverse=size(find(abs(DNAx-flip(DNAy))>0),2);
else
    error('����DNA���Ȳ�ͬ');
end
end

