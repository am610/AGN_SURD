"""Independent audit checks. Existing research outputs are read only."""
from pathlib import Path
import sys, json, hashlib, re
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'SURD/utils'))
import surd
import run_standard_iccf as iccf
import run_knn_estimator_crosscheck as knn
P=ROOT/'agn_surd_project/processed'
report={}
cont,lines,lo,hi=iccf.load_native_light_curves()
grid=np.arange(lo,hi+1)
report['data']={'window':[lo,hi],'grid_count':len(grid),'continuum_rows':len(cont),'line_rows':len(lines),'distinct_line_dates':lines.mjd.nunique(),'duplicate_line_rows':int(lines.mjd.duplicated().sum()),'fractional_continuum_dates':int((cont.time%1!=0).sum())}
for name,t in [('continuum',cont.time),('line',lines.mjd)]:
 t=np.unique(t); gaps=np.diff(t)
 bracket=np.searchsorted(t,grid,side='right');valid=(bracket>0)&(bracket<len(t)); widths=np.zeros(len(grid));widths[valid]=t[bracket[valid]]-t[bracket[valid]-1]
 widths[np.isin(grid,t)]=0
 report['data'][name]={'max_gap':float(gaps.max()),'median_gap':float(np.median(gaps)),'grid_points_bridging_gap_over_30d':int((widths>30).sum()),'grid_points_bridging_gap_over_60d':int((widths>60).sum())}
series=knn.load_real_series()
cache=np.load(P/'round_robin_surrogate_scans.npz')
report['null_recalculation']=[]
for i,name in enumerate(cache['targets']):
 real=cache['real_scans'][i];null=cache['null_scans'][:,i,:]
 report['null_recalculation'].append({'target':str(name),'peak_lag':int(cache['lags'][real.argmax()]),'maximum':float(real.max()),'p_global':float((1+(null.max(axis=1)>=real.max()).sum())/(1+len(null)))})
rows=[]; max_balance=0
for idx,name in enumerate(knn.TARGETS):
 for lag in range(1,201):
  data=np.column_stack([series[idx,lag:],series[[j for j in range(4) if j!=idx],:-lag].T]);h=np.histogramdd(data,bins=8)[0];empty=float((h==0).mean());r,s,mi,leak=surd.surd(h.copy());joint=mi[(1,2,3)];syn=sum(s.values());balance=abs(sum(r.values())+syn-joint);max_balance=max(max_balance,balance)
  rows.append({'target':name,'lag':lag,'synergy_bits':syn,'joint_mi_bits':joint,'normalized_synergy':syn/joint,'leak':leak,'empty_fraction':empty})
fresh=pd.DataFrame(rows);fresh.to_csv(ROOT/'audit/recomputed_scans.csv',index=False)
report['fresh_scan_max_difference']={}
for i,name in enumerate(cache['targets']):
 v=fresh[fresh.target==name].normalized_synergy.to_numpy();report['fresh_scan_max_difference'][str(name)]=float(np.max(np.abs(v-cache['real_scans'][i])))
report['max_atom_additivity_error']=max_balance
report['iccf_fresh']=[]
for name,col in [('blue_wing','blue_wing_flux'),('core','core_flux'),('red_wing','red_wing_flux')]:
 c=iccf.bidirectional_iccf(cont.time.to_numpy(),cont.flux.to_numpy(),lines.mjd.to_numpy(),lines[col].to_numpy());peak,cent,r=iccf.peak_and_centroid(iccf.LAGS,c)
 saved=pd.read_csv(P/'iccf_curves.csv');saved=saved[saved.component==name].r.to_numpy()
 dedup=lines.groupby('mjd',as_index=False).mean(numeric_only=True)
 c2=iccf.bidirectional_iccf(cont.time.to_numpy(),cont.flux.to_numpy(),dedup.mjd.to_numpy(),dedup[col].to_numpy());peak2,cent2,r2=iccf.peak_and_centroid(iccf.LAGS,c2)
 report['iccf_fresh'].append({'target':name,'peak':peak,'centroid':cent,'rmax':r,'stored_curve_max_difference':float(np.max(abs(c-saved))),'duplicate_dates_mean_centroid':cent2})
# Exact probability distributions test the atom implementation independently of AGN assumptions.
cases={}
for mode in ['unique','redundancy','xor','independent']:
 p=np.zeros((2,2,2))
 for a in range(2):
  for b in range(2):
   if mode=='unique':p[a,a,b]+=.25
   elif mode=='redundancy':p[a,a,a]+=.25
   elif mode=='xor':p[a^b,a,b]+=.25
   else:
    for y in range(2):p[y,a,b]+=.125
 r,s,mi,l=surd.surd(p);cases[mode]={'U1':r[(1,)],'U2':r[(2,)],'R12':r[(1,2)],'S12':s[(1,2)],'leak':l}
report['analytic_atoms']=cases
# Compare each saved spectrum reduction with its originating spectrum.
spectral=[]
for p in sorted((ROOT/'agn_surd_project/agn_data/ngc5548_agnwatch/hb_profiles_extracted').glob('*.spc')):
 a=np.loadtxt(p);v=(a[:,0]/(4861.33*1.017175)-1)*299792.458
 d={'filename':p.name,'time':float(re.search(r'n(\d{5})',p.name)[1])-10000,'pixel_width_min':float(np.diff(a[:,0]).min()),'pixel_width_max':float(np.diff(a[:,0]).max())}
 for name,low,high in [('blue_wing',-6000,-2000),('core',-2000,2000),('red_wing',2000,6000)]:
  m=(v>=low)&(v<high);d[name+'_flux']=a[m,1].sum();d[name+'_error']=np.sqrt((a[m,2]**2).sum())
 spectral.append(d)
spec=pd.DataFrame(spectral);spec.to_csv(ROOT/'audit/spectral_reduction_check.csv',index=False)
report['spectra']={'count':len(spec),'pixel_width_min':spec.pixel_width_min.min(),'pixel_width_max':spec.pixel_width_max.max()}
cols=[c for c in lines if c.endswith('_flux') or c.endswith('_error')]
a=spec.sort_values(['time']+cols)[cols].to_numpy();b=lines.sort_values(['mjd']+cols)[cols].to_numpy();report['spectra']['saved_flux_error_max_difference']=float(np.max(np.abs(a-b)))
report['source_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'SURD/utils/surd.py',ROOT/'run_standard_iccf.py',ROOT/'run_knn_estimator_crosscheck.py',P/'ngc5548_hb_velocity_bins.csv',P/'round_robin_surrogate_scans.npz',ROOT/'overleaf_draft/final_draft.tex']}
(ROOT/'audit/verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
