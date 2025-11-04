function Ydot_current = Ydot(Y_current,Y_previous,Ydot_previous,gt,dt)

Ydot_current = (Y_current - Y_previous)./(gt*dt) - ((1-gt)/gt).*Ydot_previous;
