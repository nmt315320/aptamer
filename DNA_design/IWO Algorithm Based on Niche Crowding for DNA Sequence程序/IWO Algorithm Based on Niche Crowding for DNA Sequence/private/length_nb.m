function LenNoBlank= length_nb(IntDNA_nb)
%������ĸ��{A C G T -}->{1 2 3 4 5}���ַ����з�-���ַ���
LenNoBlank=sum(IntDNA_nb~=5);
end

