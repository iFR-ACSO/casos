function [] = plotBoxCon(dims,x_ups,x_lows)
% box constraints are scaled
line([x_lows(dims(1)) ,x_ups(dims(1))] ,[x_ups(dims(2))  ,x_ups(dims(2))],'Color','black','LineStyle','--')
line([x_lows(dims(1)) ,x_ups(dims(1))] ,[x_lows(dims(2)) ,x_lows(dims(2))],'Color','black','LineStyle','--')
line([x_lows(dims(1)) ,x_lows(dims(1))],[x_lows(dims(2)) ,x_ups(dims(2))],'Color','black','LineStyle','--')
line([x_ups(dims(1))  ,x_ups(dims(1))] ,[x_lows(dims(2)) ,x_ups(dims(2))],'Color','black','LineStyle','--')
end