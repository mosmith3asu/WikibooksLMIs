clear all
clc
%https://pdf.sciencedirectassets.com/271426/1-s2.0-S0005109806X03242/1-s2.0-S0005109806001488/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEOr%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQDzzzjAiWDsPlqM2Wmvu3gY67QBeJQgrarg1HMdm0N1eQIgAqIP3gIKRpzhbMTSoS059yPeCENgJ0T9eIpCjrVzausqvQMIg%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARADGgwwNTkwMDM1NDY4NjUiDJTMk3KvGmwxRzQ6dyqRA0naya7CDAEoerzRCBmEiEYTpt6q3t%2FxhONTUg5eXv8G28kEjF%2F6GHT9iq%2FypDZ6iGKc%2BvM2H%2Fd84VJsca5kXv1BYYoivmLwgJuvxExHQH1iDIvpu8p8lwnXKTx6g3LKIUu7e0IZDhuTjPQjh2gYT%2F0wXHPDgTTWNA8upc8mTYm9FegG5XcKdQ%2BejUbp1OprbnBGFuWynksgaNFWQAzJJhpUfKdEHpXVNmbVf9b5oeAeefps3DTgZivlevhCylyC%2FA2AJTsLKkPkP9sXTZPW6QPLvHpUeeFjqhMHaZWZAnZfmBUF%2Fs1JdbqvDYMGFAKmRn%2BBnOt00LMRj1tu1qr0tlOhFhKlXC4qPszGkndpOf%2BVmOgB%2FdVdIW3wUM%2FAYL%2FCQdlhETbc1hTwlt5xHv1kbmP0X%2FWw5z%2Bmd5JT%2FfjtJh4Zz39FwzWVBC%2FddT2dmwEPfINRzo58xyjaOrFgG21qs6Bwy4GvwukNK%2BcPnR1LUqY52c9NYPc1P82NIaJ51aaCBI8nUWho1awB8Gor%2Bbxj3Oi4MP%2Bly%2F4FOusBJDc74is1E%2BAYsr2xaVw3pqp5nEfyB%2FjudFJwUPByJRiZemPaEWhTi6sRl3OuKfpyuAwSF8paqb7%2F3DogrtDQrX90hvAdlBk2SF9n7r3a5PXQFb8eRQ2h%2BWORNI5bXDkTVDEyiQK8f1Jlm6uN7E5yB5aPJhU3MTFvk4kFuqAs2uW9NNJtAomO0g7aR%2FiaroYhw1%2BH60dXtpEZmKjuj%2FpZqWPZ6dWcDBk%2BkZGgl9CbEdpcRRgry%2BuqH2GMErhl11nCoin9Yp53%2BvEa%2BhtlMtzyf%2F42adQi0uOOa9alNMXuNt8dMK%2Fov3%2FPNR7oVg%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20201211T031122Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY7BUA3T5C%2F20201211%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=2800ae916fc4c4637a98e033be3f9f81b0449342bc957a8548e4bde8515f289d&hash=f3dbb14955b50542d84c797d1451f260ac9e41b0acf30f4400c5e5bd6d61886e&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0005109806001488&tid=spdf-4737e519-5ae9-4c39-bdb7-56d37b29661c&sid=73c16225879aa2456e88bd84b6ebe22bd225gxrqa&type=client
options = sdpsettings('solver','mosek','verbose',0);
pi=[0.01 0 0.01];
d_hat=1;
y_hat=100;
wmax=10;
eps=0.00001;
for i =1:2
    Ai=[0 0; 0 -pi(i)];
    bi=[0; -d_hat];
    Bi = [1 -1; 0 1];
    Bwi=[0; -1];
    Ci=[0 1];
    %Form 11
    Si=sdpvar(2);
    Q=sdpvar(2);
    eta=sdpvar(1);

    MAT = [
        Si+Si'+eta*Bwi*Bwi' Q*Ci';
        Ci*Q eye(1)];
    
    F=[];
    F=[F;MAT<=0];
    
    F=[F;Q>=eps*eye(size(Q))];
    F=[F;eta>=eps];
    
    %F=[F;Si>=eps*eye(size(Si))];
    optimize(F,-eta,options);
    
    Q=value(Q);
    Si=value(Si);
    gamma=sqrt(value(eta)) %eta = gamma^2
    %Si=Ai*Q+Bi*Yi
    Yi=(Si+Ai*Q)*inv(Bi);
    
    K= Yi*Q^-1
end

