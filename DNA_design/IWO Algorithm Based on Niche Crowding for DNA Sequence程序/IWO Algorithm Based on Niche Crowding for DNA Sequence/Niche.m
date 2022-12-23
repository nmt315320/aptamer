function lp = Niche(popsize,chromlength,L,penalty,lp)
for i = 1:popsize-1
    for j = i+1:popsize
        d(i,j) = hamming(chromlength,lp(i,1:chromlength),lp(j,1:chromlength));%�������ĺ�������
        if d(i,j) <= L %�Ƚϸ���͸������Ӧ�ȴ�С������������Ӧ�Ƚϵ͵ĸ��崦�Է�����
            if (lp(i,chromlength+1) <= lp(j,chromlength+1))
                lp(j,chromlength+1) = penalty;
            else
                lp(i,chromlength+1) = penalty;
            end
        else
            continue;
        end
    end
end
