function LenNoBlank= length_nb(IntDNA_nb)
%返回字母表{A C G T -}->{1 2 3 4 5}上字符串中非-的字符数
LenNoBlank=sum(IntDNA_nb~=5);
end

