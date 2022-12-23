function distHammingReverse = dHR( DNAx,DNAy )
%计算DNAx与DNAy的逆序列的汉明距离,自变量对称函数
xLen=size(DNAx,2);
yLen=size(DNAy,2);
if xLen==yLen
    distHammingReverse=size(find(abs(DNAx-flip(DNAy))>0),2);
else
    error('两个DNA长度不同');
end
end

