function BoolValue = bp( IntDNAa,IntDNAb )
%intDNAa与intDNAb是否互补，Y返回1,N返回0;IN：整数编码DNA向量{ACGT}->{2031};OUT：整数值{0,1}
l=size(IntDNAa,2);
m=size(IntDNAb,2);
if l==m
    temp=3-IntDNAa;
    if sum(temp==IntDNAb)==l
        BoolValue=1;%'互补'
    else
        BoolValue=0;%'不互补'
    end
else
    BoolValue=0;
    error('DNA长度不同');
end

