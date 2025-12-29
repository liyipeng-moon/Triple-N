function add_onset_area(ratio, onset_time)

    yl=get(gca,'YLim');
    y2 = yl(1)+(yl(2)-yl(1))/ratio;
    y1 = yl(1);
    patch([0,onset_time,onset_time,0],[y1,y1,y2,y2],'k','FaceAlpha',0.15,'EdgeAlpha',0)

end