function BoolValue = bp( IntDNAa,IntDNAb )
%intDNAa��intDNAb�Ƿ񻥲���Y����1,N����0;IN����������DNA����{ACGT}->{2031};OUT������ֵ{0,1}
l=size(IntDNAa,2);
m=size(IntDNAb,2);
if l==m
    temp=3-IntDNAa;
    if sum(temp==IntDNAb)==l
        BoolValue=1;%'����'
    else
        BoolValue=0;%'������'
    end
else
    BoolValue=0;
    error('DNA���Ȳ�ͬ');
end

