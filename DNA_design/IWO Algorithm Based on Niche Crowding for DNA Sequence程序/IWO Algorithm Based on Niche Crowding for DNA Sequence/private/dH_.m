function CommonBases = dH_( DNAx,DNAy )
%计算DNAx与DNAy的相似碱基数,自变量对称函数
xLen=size(DNAx,2);
yLen=size(DNAy,2);
if xLen==yLen
CommonBases=size(find(abs(DNAx-DNAy)==0),2);
else
    error('两个DNA长度不同');
end
end


