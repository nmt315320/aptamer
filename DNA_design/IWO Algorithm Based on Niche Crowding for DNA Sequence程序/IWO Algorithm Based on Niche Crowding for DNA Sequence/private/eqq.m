function BoolValue= eqq(IntDNAa,IntDNAb)
%intDNAa��intDNAb�Ƿ���ͬ��Y����1,N����0;IN����������DNA����{ACGT}->{2031};OUT������ֵ{0,1}
if size(IntDNAa,2)==size(IntDNAb,2)
    if sum(IntDNAa==IntDNAb)==size(IntDNAa,2)
        BoolValue=1;% '��ͬ'
    else
        BoolValue=0;% '��ͬ'
    end
else
    BoolValue=0;
    error('DNA���Ȳ�ͬ');    
end
end

