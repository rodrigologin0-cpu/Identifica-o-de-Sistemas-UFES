"""MOSTURACAO 2026.3 — curvas completas de receitas 01 e 02 do artigo IFES (2018).

Compara: PID fixo, PID por faixas, Fuzzy-PID Gain Scheduling e Fuzzy puro.
Modelos GP1..GP4 de 2018 interpolados conforme temperatura medida, como
HIPOTESE LPV. Unidade de MV assumida: pontos percentuais de 0-100%.
NAO e validacao fisica da panela real; duracoes de patamar sao HIPOTESES.

Requisitos: numpy, matplotlib. numba e opcional.
Execucao: python mosturacao_2026_curvas_completas.py

Pode reotimizar sobre as DUAS RECEITAS inteiras; os parametros de
2026.2 sao apenas pontos de partida, nao conclusao final.
"""
from pathlib import Path
import csv, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    from numba import njit
except ImportError:
    def njit(f): return f

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'mosturacao_2026_definitivo_resultados'
OUT.mkdir(exist_ok=True,parents=True)
PARAM_FILE = OUT / 'parametros_definitivos.json'
# Mesmos ajustes publicados no pacote comparativo 2026.2, embutidos para execucao independente.
FALLBACK_PARAMS = [
    [7.032974589986552, 526.3905843177542, 15.343351746438946],
    [0.7684252833699396, 1.3469023512619822, 0.6512309662059637],
    [0.9183393191084079, 1.4317599759399116, 0.6396130453516156, 0.8175470113129372, 0.7575335905456689, 0.5092177994730611],
    [0.24667151331596837, 2.2181656665493024, 3.803819746131396, 0.14228038821946576],
]
TAMB = 26.0
DT = 1.0
MAX_MV = 100.0
BAND_C = 0.2
NUMERIC_OVERSHOOT_TOL_C = 0.02  # tolerancia numerica de relatorio, NAO banda de qualidade
MIN_CONTINUOUS_HOLD = 600       # 10min na banda, hipotese de comparacao
D_FILTER_S = 10.0
MODEL_K = np.array([1.129, 1.016, 1.011, .981],dtype=float)
MODEL_TAU = np.array([366., 292., 216., 178.],dtype=float)
MODEL_DELAY = np.array([28., 19., 14., 8.7],dtype=float)
TEMP_CENTERS = np.array([42., 52., 63., 70.5],dtype=float)  # hipotese, NAO identificacao de 2018
CHR = np.array([[6.94,366.,14.],[8.36,292.,9.5],[9.15,216.,7.],[12.51,178.,4.35]],dtype=float)
CONTROL_LABELS = ['PID fixo','PID Gain Scheduling por faixas','Fuzzy-PID Gain Scheduling','Fuzzy puro incremental']
RECIPES = {
    'Receita 01': {
        'temps':[62.,70.5,78.],
        'labels':['Beta-amilase','Alfa-amilase','Inativacao / final'],
        'durations':[1800,1200,1200],
    },
    'Receita 02': {
        'temps':[44.,52.,63.,70.5,78.],
        'labels':['Beta-glucanase','Parada proteica','Beta-amilase','Alfa-amilase','Inativacao / final'],
        'durations':[1800,1200,1200,1200,1200],
    }
}

@njit
def clip(v,lo,hi):
    return min(max(v,lo),hi)

@njit
def linear_piecewise(x,ordinates):
    if x<=TEMP_CENTERS[0]:return ordinates[0]
    if x>=TEMP_CENTERS[-1]:return ordinates[-1]
    for j in range(3):
        if x<TEMP_CENTERS[j+1]:
            a=(x-TEMP_CENTERS[j])/(TEMP_CENTERS[j+1]-TEMP_CENTERS[j])
            return (1.-a)*ordinates[j]+a*ordinates[j+1]
    return ordinates[-1]

@njit
def region_index(t):
    if t<47:return 0
    if t<57.5:return 1
    if t<66.75:return 2
    return 3

@njit
def memberships4(temp):
    w=np.zeros(4)
    if temp<=TEMP_CENTERS[0]:w[0]=1.;return w
    if temp>=TEMP_CENTERS[3]:w[3]=1.;return w
    for j in range(3):
        if TEMP_CENTERS[j]<=temp<TEMP_CENTERS[j+1]:
            f=(temp-TEMP_CENTERS[j])/(TEMP_CENTERS[j+1]-TEMP_CENTERS[j])
            w[j]=1-f;w[j+1]=f
            return w
    w[3]=1.
    return w

@njit
def memberships3_error(e):
    w=np.zeros(3)
    if e<=0.:w[0]=1.
    elif e<.8:w[0]=1.-e/.8;w[1]=e/.8
    elif e<3.:w[1]=(3.-e)/(3.-.8);w[2]=(e-.8)/(3.-.8)
    else:w[2]=1.
    return w

@njit
def fuzzy_gain(temp,e,p):
    tw=memberships4(temp);ew=memberships3_error(e)
    kp=0.;ki=0.;kd=0.;total=0.
    for j in range(4):
        for k in range(3):
            w=min(tw[j],ew[k])
            if w<=0:continue
            scale=p[3] if k==0 else (1.0 if k==1 else p[4])
            kt=CHR[j,0]*p[0]*scale
            it=CHR[j,1]*p[1]
            td=CHR[j,2]*p[2]*(p[5] if k==0 else 1.)
            kp+=w*kt;ki+=w*kt/it;kd+=w*kt*td;total+=w
    if total<1e-10:return 8.,8./292.,70.
    return kp/total,ki/total,kd/total

@njit
def triangle5(x):
    w=np.zeros(5)
    if x<=-2:w[0]=1.;return w
    if x>=2:w[4]=1.;return w
    j=int(math.floor(x+2.))
    f=x-(-2.+j)
    w[j]=1.-f;w[j+1]=f
    return w

@njit
def fuzzy_increment(error,vel,p):
    e=clip(error/p[2],-2.,2.)
    d=clip(vel/p[3],-2.,2.)
    ew=triangle5(e);dw=triangle5(d)
    num=0.;den=0.
    for i in range(5):
        for j in range(5):
            w=min(ew[i],dw[j]);
            if w<=0.:continue
            consequence=.65*(i-2)-p[1]*(j-2)
            num+=w*consequence;den+=w
    if den<=0:return 0.
    return clip(p[0]*num/den,-4.,4.)

@njit
def simulate(which,p,sp_array,kfac=1.,tfac=1.,extra_delay=0.):
    n=len(sp_array)
    y=np.zeros(n); u=np.zeros(n)
    y[0]=TAMB
    gains=np.zeros((n,3))
    i_term=0.;vfilter=0.
    for k in range(n):
        temp=y[k]
        err=sp_array[k]-temp
        if k>0:
            inst=(y[k]-y[k-1])/DT
            alpha=math.exp(-DT/D_FILTER_S)
            vfilter=alpha*vfilter+(1.-alpha)*inst
        if which==0:
            kp=p[0];ki=kp/p[1];kd=kp*p[2]
        elif which==1:
            idx=region_index(temp)
            kp=CHR[idx,0]*p[0]
            ki=kp/(CHR[idx,1]*p[1]);kd=kp*CHR[idx,2]*p[2]
        elif which==2:
            kp,ki,kd=fuzzy_gain(temp,abs(err),p)
        else:
            kp=0.;ki=0.;kd=0.
        gains[k,0]=kp;gains[k,1]=ki;gains[k,2]=kd
        if which!=3:
            prop=kp*err
            tentative_i=i_term+ki*err*DT
            trial=prop+tentative_i-kd*vfilter
            if not ((trial>MAX_MV and err>0.) or (trial<0. and err<0.)):
                i_term=tentative_i
            u[k]=clip(prop+i_term-kd*vfilter,0.,MAX_MV)
        else:
            prev=u[k-1] if k>0 else 0.
            u[k]=clip(prev+fuzzy_increment(err,vfilter,p),0.,MAX_MV)
        if k==n-1:break
        # Planta LPV com parametros interpolados por temperatura corrente.
        # EQUILIBRIO: 26°C + K(T) * MV. Hipotese, nao calibracao fisica.
        gain=linear_piecewise(temp,MODEL_K)*kfac
        tau=linear_piecewise(temp,MODEL_TAU)*tfac
        theta=linear_piecewise(temp,MODEL_DELAY)+extra_delay
        pos=k-theta/DT
        ix=int(math.floor(pos))
        frac=pos-ix
        delayed=0.
        if ix>=0 and ix+1<n:
            delayed=(1-frac)*u[ix]+frac*u[ix+1]
        elif ix>=0:
            delayed=u[ix]
        a=math.exp(-DT/tau)
        y[k+1]=TAMB+(y[k]-TAMB)*a+gain*delayed*(1.-a)
    return y,u,gains


def setpoint_schedule(r):
    durations=r['durations'];temps=r['temps']
    starts=np.concatenate(([0],np.cumsum(durations)))
    arr=np.zeros(int(starts[-1])+1)
    for i,t in enumerate(temps):
        arr[int(starts[i]):int(starts[i+1])]=t
    arr[-1]=temps[-1]
    return arr,starts


def metrics(y,u,sp,r,which,recipe,case):
    rows=[]
    durations=r['durations']
    starts=np.concatenate(([0],np.cumsum(durations)))
    for i,tar in enumerate(r['temps']):
        beg=int(starts[i]);end=int(starts[i+1]);seg=y[beg:end]
        mv=u[beg:end]
        diff=seg-tar
        mp=max(0.,float(np.max(diff)))
        ov=diff>NUMERIC_OVERSHOOT_TOL_C
        overs_time=int(ov.sum())
        area=float(np.maximum(diff,0.).sum()*DT)
        # Tempo para primeiro bloco CONTINUO de 600 segundos dentro da banda.
        inband=np.abs(diff)<=BAND_C
        longest=0;counter=0;first_hold=None
        for ix,b in enumerate(inband):
            counter=counter+1 if b else 0
            longest=max(longest,counter)
            if counter>=MIN_CONTINUOUS_HOLD and first_hold is None:
                first_hold=ix+1-MIN_CONTINUOUS_HOLD
        # Ts real (critério sem saída posterior da banda ate proximo SP)
        last_bad=np.where(~inband)[0]
        ts_full = (int(last_bad[-1]+1) if len(last_bad) else 0)
        if ts_full>=len(seg):ts_full=None
        if ts_full is not None and (len(seg)-ts_full)<MIN_CONTINUOUS_HOLD:ts_full=None
        rows.append({
            'receita':recipe,'cenario':case,'controlador':CONTROL_LABELS[which],
            'patamar':i+1,'fase':r['labels'][i], 'SP_C':tar,
            'inicio_s':beg,'fim_s':end,'duracao_s':end-beg,
            'T_inicial_C':round(float(seg[0]),4),'T_final_C':round(float(seg[-1]),4),
            'Ts_s':ts_full if ts_full is not None else '',
            'primeiro_bloco_600s_s':first_hold if first_hold is not None else '',
            'tempo_max_continuo_banda_s':longest,
            'sustentou_10min':bool(longest>=MIN_CONTINUOUS_HOLD),
            'sobressinal_C':round(mp,5),
            'sem_sobressinal_0_02C':bool(mp<=NUMERIC_OVERSHOOT_TOL_C),
            'tempo_acima_SP_s':overs_time,
            'area_excesso_C_s':round(area,3),
            'tempo_na_banda_s':int(inband.sum()),
            'MV_media_pct_assumida':round(float(mv.mean()),3),
        })
    return rows



# --- ETAPA FINAL: SINTONIA COMPARATIVA EM CURVAS COMPLETAS ---
import argparse
from scipy.optimize import differential_evolution

# Cenarios de treino e validacao sao separados. Centros das regioes,
# normalizacao K e duracoes continuam HIPOTESES declaradas.
TRAIN_CASES=[('nominal',1.0,1.0,0.0),
             ('treino_mais_lento',0.94,1.08,3.0),
             ('treino_mais_rapido',1.07,0.93,1.0),
             ('treino_ganho_alto',1.12,0.90,5.0),
             ('treino_ganho_baixo',0.88,1.15,7.0)]
TEST_CASES=[('teste_ganho_alto',1.15,0.85,6.0),
            ('teste_lento',0.85,1.20,9.0),
            ('teste_atraso',1.05,1.0,12.0),
            ('teste_k_baixo_tau_rapido',0.90,0.87,8.0)]
# Restricao numerica, a resolucao temporal e de 1 segundo. O controle real
# exige margens adicionais e nao e validado pelo criterio abaixo.
MAX_ACCEPTED_OVERSHOOT=0.02

@njit
def quality_cost(y,targets,durations):
    start=0
    n_over=0
    max_over_sum=0.
    n_failed=0
    sum_ts=0.
    sum_remaining_gap=0.
    for j in range(len(durations)):
        dur=int(durations[j])
        stop=start+dur
        tar=targets[j]
        maxover=0.
        lastbad=-1
        for k in range(start,stop):
            delta=y[k]-tar
            if delta>maxover: maxover=delta
            if abs(delta)>BAND_C: lastbad=k-start
        if maxover>MAX_ACCEPTED_OVERSHOOT:
            n_over+=1
            max_over_sum+=maxover-MAX_ACCEPTED_OVERSHOOT
        ts=lastbad+1
        if dur-ts<MIN_CONTINUOUS_HOLD:
            n_failed+=1
            sum_remaining_gap+=max(0.0,MIN_CONTINUOUS_HOLD-(dur-ts))
            ts=dur+MIN_CONTINUOUS_HOLD
        sum_ts+=ts
        start=stop
    return n_over,max_over_sum,n_failed,sum_ts,sum_remaining_gap

SP1,_=setpoint_schedule(RECIPES['Receita 01'])
SP2,_=setpoint_schedule(RECIPES['Receita 02'])
TARGET1=np.array(RECIPES['Receita 01']['temps'],dtype=float)
TARGET2=np.array(RECIPES['Receita 02']['temps'],dtype=float)
DURS1=np.array(RECIPES['Receita 01']['durations'],dtype=np.int64)
DURS2=np.array(RECIPES['Receita 02']['durations'],dtype=np.int64)

@njit
def objective(p,which):
    overs=0; severity=0.; missing=0; time_sum=0.; lost=0.
    for ridx in range(2):
        if ridx==0: sp=SP1;temps=TARGET1;durations=DURS1
        else: sp=SP2;temps=TARGET2;durations=DURS2
        for k in range(5):
            if k==0: kf,tf,delay=1.0,1.0,0.0
            elif k==1: kf,tf,delay=.94,1.08,3.
            elif k==2: kf,tf,delay=1.07,.93,1.
            elif k==3: kf,tf,delay=1.12,.90,5.
            else: kf,tf,delay=.88,1.15,7.
            y,_,_=simulate(which,p,sp,kf,tf,delay)
            ov,sev,miss,ts,lg=quality_cost(y,temps,durations)
            overs+=ov;severity+=sev;missing+=miss;time_sum+=ts;lost+=lg
    # Ordem de prioridade explícita: 1) sobressinal zero, 2) manter patamar,
    # 3) reduzir Ts. Punicao continua dentro de cada nivel para orientar busca.
    if overs>0:
        return 2.e7 + 7.e5*overs + 7.e6*severity + 4.e4*missing + .001*time_sum
    if missing>0:
        return 2.e6 + 5.e4*missing + 10.*lost + .01*time_sum
    return time_sum / 40.0

BOUNDS=[[(1.0,18.0),(50.,1100.),(0.,70.)],
        [(.3,2.5),(.3,2.8),(.0,4.)],
        [(.3,2.5),(.3,2.8),(.1,4.),(.15,1.4),(.55,2.5),(.3,3.)],
        [(.03,3.3),(.1,3.0),(.15,6.),(.008,.25)]]


def optimize(maxiter=25,popsize=5):
    params=[];reports=[]
    for which,name in enumerate(CONTROL_LABELS):
        print(f'Sintonia [{which+1}/4] {name}',flush=True)
        # Com as mesmas janelas, ambientes e numero de geracoes. Maior numero
        # de parametros = mais graus de liberdade, nao conclusao universal.
        res=differential_evolution(
            lambda z:objective(np.asarray(z,dtype=float),which),
            BOUNDS[which],maxiter=maxiter,popsize=popsize,seed=2026+which,
            polish=False,tol=0.006,init='latinhypercube',workers=1,
            updating='immediate',disp=False)
        params.append(res.x.tolist())
        reports.append({'controlador':name,'custo_treino':float(res.fun),
                        'iteracoes':res.nit,'avaliacoes':res.nfev,'sucesso_otimizador':bool(res.success)})
        print('  custo:',round(res.fun,2),'param:',np.round(res.x,5),'nfev:',res.nfev,flush=True)
    package={'versao':'2026.4 - curvas completas; ajuste igualitário',
             'parametros':params,'relatorio_otimizacao':reports,
             'semente':2026,'maxiter':maxiter,'popsize':popsize,
             'criterio':'Prioridade 1 limite de sobressinal 0.02 na escala hipotetica; prioridade 2 600s consecutivos em +/-0.2; prioridade 3 menor Ts medio; 5 cenarios de treino separados dos 4 de validacao',
             'limites':'0..100% assumidos, GP1..GP4 com centros de temperatura hipoteticos; simulacao apenas'}
    PARAM_FILE.write_text(json.dumps(package,indent=2,ensure_ascii=False),encoding='utf-8')
    return [np.array(x,dtype=float) for x in params]


def make_final_plots(recipe,sp,starts,sims,scenario):
    slug='receita_01' if '01' in recipe else 'receita_02'
    suffix='nominal' if scenario=='nominal' else scenario
    tm=np.arange(len(sp))/60
    fig,axes=plt.subplots(2,1,figsize=(15,8.6),sharex=True,layout='constrained',height_ratios=[2.5,1])
    axes[0].step(tm,sp,label='Setpoint (receita do artigo; durações hipotéticas)',where='post',color='#333333',ls='--',lw=2.2)
    colors=['#555555','#C98210','#087E8B','#7654A3']
    for idx,label in enumerate(CONTROL_LABELS):
        y,u,_=sims[idx]
        axes[0].plot(tm,y,color=colors[idx],label=label,lw=1.6)
        axes[1].plot(tm,u,color=colors[idx],label=label,lw=1.3)
    axes[0].set(ylabel='Temperatura (°C, modelo hipotético)',title=f'Mosturação 2026.4 | {recipe} | {scenario}')
    axes[0].legend(fontsize=9,loc='lower right')
    axes[1].set(ylabel='MV (% assumido)',xlabel='Tempo (min)',ylim=(-3,103))
    for ax in axes:
        for s in starts[1:-1]:ax.axvline(s/60,c='#999999',ls=':',alpha=.32)
        ax.grid(alpha=.16)
    fig.savefig(OUT/f'{slug}_{suffix}_curva.png',dpi=155)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(15,5.3),layout='constrained')
    for idx,label in enumerate(CONTROL_LABELS):
        ax.plot(tm,sims[idx][0]-sp,color=colors[idx],lw=1.5,label=label)
    for yy,ls,col in [(0,'--','#111111'),(BAND_C,':','#888888'),(-BAND_C,':','#888888')]:ax.axhline(yy,ls=ls,c=col,lw=1.)
    ax.set(xlabel='Tempo (min)',ylabel='PV – SP (°C, escala hipotética)',title=f'Erro e ultrapassagem | {recipe} | {scenario}')
    ax.legend(fontsize=9);ax.grid(alpha=.18)
    fig.savefig(OUT/f'{slug}_{suffix}_erro.png',dpi=155)
    plt.close(fig)


def save_csv(filename,rows):
    if not rows:return
    with (OUT/filename).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter=';')
        w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description='Mosturação 2026.4; comparativo com curvas completas')
    parser.add_argument('--otimizar',action='store_true',help='Refazer sintonia dos quatro controladores')
    parser.add_argument('--maxiter',type=int,default=25)
    parser.add_argument('--popsize',type=int,default=5)
    opt=parser.parse_args()
    if opt.otimizar or not PARAM_FILE.exists():
        params=optimize(maxiter=opt.maxiter,popsize=opt.popsize)
    else:
        data=json.loads(PARAM_FILE.read_text(encoding='utf-8'))
        params=[np.array(z,dtype=float) for z in data['parametros']]
    rows=[];samples=[]
    for recipe,r in RECIPES.items():
        sp,starts=setpoint_schedule(r)
        for group, cases in [('treino',TRAIN_CASES),('validacao',TEST_CASES)]:
            for sc,kf,tf,delay in cases:
                curves={}
                for idx in range(4):
                    y,u,g=simulate(idx,params[idx],sp,kf,tf,delay)
                    curves[idx]=(y,u,g)
                    rr=metrics(y,u,sp,r,idx,recipe,sc)
                    for item in rr:
                        item['grupo']=group
                        item['fator_K']=kf;item['fator_tau']=tf;item['atraso_extra_s']=delay
                        item['atingiu_sem_ultrapassagem']=bool(item['sustentou_10min'] and item['sobressinal_C']<=MAX_ACCEPTED_OVERSHOOT)
                    rows.extend(rr)
                    if sc=='nominal':
                        for k in range(len(sp)):
                            samples.append({'receita':recipe,'controlador':CONTROL_LABELS[idx],
                                'segundo':k,'SP_C':sp[k],'PV_C':round(float(y[k]),5),
                                'MV_pct_assumida':round(float(u[k]),5),
                                'Kp':round(float(g[k,0]),5),'Ki':round(float(g[k,1]),7),'Kd':round(float(g[k,2]),5)})
                if sc=='nominal' or sc in ['teste_ganho_alto','teste_lento']:
                    make_final_plots(recipe,sp,starts,curves,sc)
    save_csv('metricas_por_patamar.csv',rows)
    save_csv('series_temporais_nominais.csv',samples)
    summary=[]
    for group,cases in [('treino',TRAIN_CASES),('validacao',TEST_CASES)]:
        for idx,lab in enumerate(CONTROL_LABELS):
            rr=[x for x in rows if x['grupo']==group and x['controlador']==lab]
            ts=[x['Ts_s'] for x in rr if x['Ts_s']!='']
            valid=[x for x in rr if x['atingiu_sem_ultrapassagem']]
            summary.append({'grupo':group,'controlador':lab,'patamares':len(rr),
                'qualidade_10min_sem_sobressinal':len(valid),
                'cumpriu_10min':sum(x['sustentou_10min'] for x in rr),
                'patamares_sem_sobressinal':sum(x['sobressinal_C']<=MAX_ACCEPTED_OVERSHOOT for x in rr),
                'Ts_mediano_s_somente_acomodados':round(float(np.median(ts)),2) if ts else '',
                'Ts_amostras':len(ts),'pior_sobressinal_C':max(x['sobressinal_C'] for x in rr),
                'area_excesso_Cs':round(sum(x['area_excesso_C_s'] for x in rr),2)})
    save_csv('resumo_comparativo.csv',summary)
    by_recipe=[]
    for recipe in RECIPES:
        for group in ('treino','validacao'):
            for lab in CONTROL_LABELS:
                rr=[x for x in rows if x['receita']==recipe and x['grupo']==group and x['controlador']==lab]
                ts=[x['Ts_s'] for x in rr if x['Ts_s']!='']
                by_recipe.append({'receita':recipe,'grupo':group,'controlador':lab,
                    'patamares_total':len(rr),
                    'patamares_aprovados':sum(x['atingiu_sem_ultrapassagem'] for x in rr),
                    'Ts_mediano_entre_acomodados_s':round(float(np.median(ts)),2) if ts else '',
                    'sobressinal_maximo_C':max(x['sobressinal_C'] for x in rr)})
    save_csv('resumo_por_receita.csv',by_recipe)
    names=['PID fixo','PID por faixas','Fuzzy-PID','Fuzzy puro']
    for grp,filename in [('treino','sintese_tempo_assentamento.png'),('validacao','sintese_validacao.png')]:
        sub=[q for q in summary if q['grupo']==grp]
        fig,ax=plt.subplots(figsize=(11.4,5.7),layout='constrained')
        xp=np.arange(4)
        if grp=='treino':
            vals=[v['Ts_mediano_s_somente_acomodados'] for v in sub]
            bars=ax.bar(xp,vals,color=['#666666','#C98210','#087E8B','#7654A3'],width=.64)
            for i,b in enumerate(bars):
                ax.text(b.get_x()+b.get_width()/2,b.get_height()+8,
                  f"{vals[i]:g} s | {sub[i]['qualidade_10min_sem_sobressinal']}/40",ha='center',fontsize=10)
            ax.set(ylabel='Ts mediano (s) — somente patamares acomodados',title='Curvas completas | Resultados da otimização, 40 patamares por controlador')
            ax.set_ylim(0,max(vals)*1.25)
        else:
            vals=[v['qualidade_10min_sem_sobressinal'] for v in sub]
            bars=ax.bar(xp,vals,color=['#666666','#C98210','#087E8B','#7654A3'],width=.64)
            for i,b in enumerate(bars):
                ax.text(b.get_x()+b.get_width()/2,b.get_height()+.5,
                  f"{vals[i]}/32 | máx {sub[i]['pior_sobressinal_C']:.3f}",ha='center',fontsize=10)
            ax.set(ylabel='Patamares aprovados em 32 testes',title='Validação em cenários não utilizados no ajuste (sobressinal máximo anotado)')
            ax.set_ylim(0,35)
        ax.set_xticks(xp,names); ax.grid(axis='y',alpha=.16)
        fig.savefig(OUT/filename,dpi=165)
        plt.close(fig)
    for r in summary:print(r,flush=True)
    print('Arquivos gerados em',OUT,flush=True)

if __name__=='__main__':main()