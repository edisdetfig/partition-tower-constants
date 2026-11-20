\\ The l-adic Hecke systems on newforms of weight k on Gamma0(6) with U2 = e2*2^((k-2)/2), U3 = e3*3^((k-2)/2).
Bsys(l,k,e2,e3,PREC)=
{
  my(mf=mfinit([6,k],0), d=mfdim(mf));
  if(d==0, return([0,0,[]]));
  my(T2=mfheckemat(mf,2), T3=mfheckemat(mf,3));
  my(B=matintersect(matker(T2-e2*2^((k-2)/2)*matid(d)), matker(T3-e3*3^((k-2)/2)*matid(d))));
  my(r=#B); if(r==0, return([d,0,[]]));
  my(Q=setminus([5,7,11,13],[l]));
  my(M=vector(#Q,i,matinverseimage(B,mfheckemat(mf,Q[i])*B)), Ml=matinverseimage(B,mfheckemat(mf,l)*B));
  my(L=sum(i=1,#Q,(2*i+1)*M[i]), f=charpoly(L), rts=polrootspadic(f,l,PREC), out=List());
  for(j=1,#rts,
    my(lam=rts[j], g=f\('x-lam), G=subst(g,'x,L), best=1, bv=oo);
    for(c=1,r, my(v=oo); for(i=1,r, if(G[i,c]!=0, v=min(v,valuation(G[i,c],l)))); if(v<bv, bv=v; best=c));
    my(vec=G[,best], i0=1, vi=oo);
    for(i=1,r, if(vec[i]!=0 && valuation(vec[i],l)<vi, vi=valuation(vec[i],l); i0=i));
    listput(out,[vector(#Q,i,(M[i]*vec)[i0]/vec[i0]), (Ml*vec)[i0]/vec[i0]]));
  [d,r,Q,Vec(out),poldegree(f)];
}
