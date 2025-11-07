hydro = readWAMIT(struct(),'rm3.out',[]);
A_mod = hydro.A(:,:,10);
A_mod(A_mod > 3e5) = NaN;

figure
imagesc(A_mod)
hold on
for num=1:11
    if num==6
        size = 2;
    else
        size = .5;
    end
    plot([.5 12.5],.5+[num num],'w','LineWidth',size)
    plot(.5+[num num],[.5 12.5],'w','LineWidth',size)
end
colorbar
