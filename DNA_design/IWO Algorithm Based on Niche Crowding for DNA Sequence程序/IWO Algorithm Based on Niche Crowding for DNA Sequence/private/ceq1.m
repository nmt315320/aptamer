function IntValue = ceq1?( IntDNAx,IntDNAy,i)
%从第i个位置开始长度为c的IntDNAx与IntDNAy的公共子串的长度;IN:IntDNAx、IntDNAy,整数编码DNA,i,起始位置;OUT:IntValue,最长公共字串长度
common=size(find(abs(IntDNAx-IntDNAy)==0),2);
abs(IntDNAx-IntDNAy)==0
common=find(abs(IntDNAx-IntDNAy)==0);
IntValue=common;
end


