function [Xt,Ut,Yt,noise,success] = simulateEnvironment(env,ctrl)
        global temp

T = round(env.time.T/env.time.dt);
Xt = zeros(4,T);
Ut = zeros(2,T-1);
Yt = zeros(4,T);
t = 1;
N0 = mvnrnd(zeros(4,1),env.noise.initial.covariance)';
Xt(:,t) = [env.initial.x;env.initial.y;env.initial.t;env.initial.c]+N0;
Xt(4,t) = min(max(Xt(4,t),-env.geom.a/2),env.geom.a/2);
Nu = zeros(2,T-1);
Nf = zeros(4,T-1);
while sqrt((Xt(1,t)-env.goal.x)^2+(Xt(2,t)-env.goal.y)^2) > 0.05 && t <= T
    switch ctrl.type
        case 1
            [Ut(:,t),dY] = cascadedP(Xt(:,t),env,ctrl);
        case 2
            [Ut(:,t),dY] = cascadedPD(Xt(:,t),Yt(:,t),env,ctrl);
    end
    Nu(:,t) = mvnrnd(zeros(2,1),env.noise.u.covariance);
    N = diag(Nu(:,t))+eye(2);
    Ut(:,t) = abs(N)*Ut(:,t);
    Nf(:,t) = mvnrnd(zeros(4,1),env.noise.f.covariance); 
    N = diag(Nf(:,t))+eye(4);
    [~,temp] = ode45(@(time,X)N*simulationModel(X,Ut(:,t),env),[0 env.time.dt],Xt(:,t));
    Xt(:,t+1) = temp(end,:)';
    Yt(:,t+1) = Yt(:,t) + dY*env.time.dt;
    if Xt(4,t+1) > env.geom.a/2
        idx = find(temp(:,4)>env.geom.a/2,1,'first');
        Xt(:,t+1) = temp(idx,:)';
        break
    end
    if Xt(4,t+1) < -env.geom.a/2
        idx = find(temp(:,4)<-env.geom.a/2,1,'first');
        Xt(:,t+1) = temp(idx,:)';
        break
    end
    if norm(Xt(1:2,t+1)) > 1
        idx = find(sqrt(sum(temp(:,1:2).^2,2))>1,1,'first');
        Xt(:,t+1) = temp(idx,:)';
        break
    end
    t = t+1;
end
success = sqrt((Xt(1,t)-env.goal.x)^2+(Xt(2,t)-env.goal.y)^2) <= 0.05;
Xt = Xt(:,1:t-1);
Yt = Yt(:,1:t-1);
Ut = Ut(:,1:t-2);
noise.N0 = N0;
noise.Nu = Nu;
noise.Nf = Nf;
end
