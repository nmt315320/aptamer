function f=fit(DNAs)
f(:,1)=Continuity(DNAs); %Continuity
f(:,2)=Hairpin(DNAs);    %Hairpin
[f(:,3),f(:,4)]=HmSm(DNAs);     
f(:,5)=f(:,1)+f(:,2)+f(:,3)+f(:,4);
end

