function IntDNAShifted= shift( IntDNA_nb,i )
%��IntDNA_nbƽ��iλ��i>0���ƣ�i<0,����;��λ��'-'(5)���
if i==0
    IntDNAShifted=IntDNA_nb;
    return;
end
l=size(IntDNA_nb,2);%DNA����
temp=zeros(1,abs(i))+5;
if (0<i) && (i<l)   %����iλ���i��'-'
    IntDNAShifted=[temp, IntDNA_nb(1:l-i)];
    return;
end
if (i<0) && (i>-l)  %����iλ���i��'-'
    IntDNAShifted=[IntDNA_nb(abs(i)+1:l),temp];
    return;
end
if abs(i)>=l %ȫ���滻Ϊ'-'
    IntDNAShifted=zeros(1,l)+5;
    return;
end
end

