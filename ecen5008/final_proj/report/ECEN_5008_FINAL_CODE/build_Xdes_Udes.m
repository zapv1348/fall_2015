function [Xdes,Udes,T] = build_Xdes_Udes(V,Gamma,Vdot,Gammadot,T)
Xdes = zeros(4,length(T));
Udes = zeros(2,length(T));
for i = 1:length(T)
    [~,Xdes(:,i),Udes(:,i),~,~] = trim_acc([V(i) Gamma(i) Vdot(i) Gammadot(i)]);
    %%fprintf('%d \n',i)
end
end