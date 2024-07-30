%% Animation of the segway

% Importing from the simulation
x = getElement(out.yout,'x').Values.Data;
th = getElement(out.yout,'phi').Values.Data;
n = length(x);
t_max = (out.SimulationMetadata.ModelInfo.StopTime);

r = wheel_radius;
L  = center_of_mass; % length of pendulum

% Position coordination
y = r;

base = [];
cart = [];
p_cir = [];
pendulum = [];

figure(1)
clf
a = int32(n / 100);

for t = 1:a:n
    delete(base);
    delete(cart);
    delete(p_cir);
    delete(pendulum);
    
    
    xlabel('X (m)');
    ylabel('Y (m)');
    title('Inverted Pendulum')
    hold on
    grid on


    base = plot([x(t)-5*r x(t)+5*r],[0 0],'k','LineWidth',2); % base line
    hold on;
    % position of pendulum
    px = x(t) + L*sin(th(t));
    py = y + L*cos(th(t));
    
    cart = viscircles([x(t) y],r,'Color','r','LineWidth',0.2); % center of Cart
    p_cir = viscircles([x(t) y],0.02,'Color',[1 0.1 0.1],'LineWidth',2.5); % Pendulum Trajectory(Circle)
    pendulum = plot([x(t) px],[y py],'b','LineWidth',2.5); % Pendulum rod
    
    xl1 = x(t)-2*r-L;
    xl2 = x(t)+2*r+L;
    yl1 = -L;
    yl2 = 2*r+2*L;
    xlim([xl1 xl2]);
    ylim([yl1 yl2]);
    axis(gca,"square");
    grid on;
    pause(0.1);
end
