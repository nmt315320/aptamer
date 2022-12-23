function dHRC= dHRC(DNAx,DNAy )
%计算DNAx与DNAy的逆补序列的汉明距离,自变量对称函数
xLen=size(DNAx,2);
yLen=size(DNAy,2);
if xLen==yLen
    dHRC=size(find(abs(DNAx-(5-flip(DNAy)))>0),2);
else
    error('两个DNA长度不同');
end
end

