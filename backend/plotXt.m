function plotXt(Xt,env)
    clf
    plot(Xt(1,:),Xt(2,:))
    hold on
    plot(env.goal.x,env.goal.y,'r*')
    hold off
end