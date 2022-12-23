function BoolValue= eqq(IntDNAa,IntDNAb)
%intDNAa与intDNAb是否相同，Y返回1,N返回0;IN：整数编码DNA向量{ACGT}->{2031};OUT：整数值{0,1}
if size(IntDNAa,2)==size(IntDNAb,2)
    if sum(IntDNAa==IntDNAb)==size(IntDNAa,2)
        BoolValue=1;% '相同'
    else
        BoolValue=0;% '不同'
    end
else
    BoolValue=0;
    error('DNA长度不同');    
end
end

