    figure;
    
    ax1 = axes;
    x1 = [1;2];
    f1 = [10 20; 0 0];
    l1 = [1 2; 0 0];
    u1 = [2 1; 0 0];
    barwitherr(cat(3,l1,u1),x1,f1);
    
    ax2 = axes;
    x2 = [2;1];
    f2 = [1 2; 0 0];
    l2 = [0.1 0.2; 0 0];
    u2 = [0.2 0.1; 0 0];
    barwitherr(cat(3,l2,u2),x2,f2);
    
    xlim(ax1,[0 3]);
    xlim(ax2,[0 3]);
    set(ax2,'Color','none','Position',get(ax1,'Position'),'YAxisLocation','right');