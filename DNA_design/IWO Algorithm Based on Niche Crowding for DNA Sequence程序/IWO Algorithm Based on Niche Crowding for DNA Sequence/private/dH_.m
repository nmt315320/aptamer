function CommonBases = dH_( DNAx,DNAy )
%����DNAx��DNAy�����Ƽ����,�Ա����Գƺ���
xLen=size(DNAx,2);
yLen=size(DNAy,2);
if xLen==yLen
CommonBases=size(find(abs(DNAx-DNAy)==0),2);
else
    error('����DNA���Ȳ�ͬ');
end
end


