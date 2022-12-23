function lp = Niche(popsize,chromlength,L,penalty,lp)
for i = 1:popsize-1
    for j = i+1:popsize
        d(i,j) = hamming(chromlength,lp(i,1:chromlength),lp(j,1:chromlength));%两个体间的海明距离
        if d(i,j) <= L %比较个体和个体的适应度大小，并对其中适应度较低的个体处以罚函数
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
