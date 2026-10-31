"""RT Plan Intelligence research UI. DICOM files are read from memory, not stored."""
from __future__ import annotations
from io import BytesIO
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from rtpi.dicom import read_rtplan
from rtpi.metrics import Aperture, calculate_beam, control_point_records

st.set_page_config(page_title="RT Plan Intelligence",page_icon="🔬",layout="wide")

def synthetic_demo():
    """Deterministic, non-clinical examples for UI testing."""
    boundaries=np.linspace(-40,40,17)
    aps=[]
    for i in range(16):
        x=np.arange(16)
        gap=7+1.8*np.sin(x*.7+i*.35)
        center=np.sin(x*.45+i*.15)*3
        aps.append(Aperture(center-gap/2,center+gap/2,boundaries,(-25.,25.),(-35.,35.)))
    cmw=np.linspace(0,1,len(aps)).tolist()
    beam=calculate_beam(1,"Синтетическая арка",aps,cmw,500).to_dict()
    beam["control_point_data"]=control_point_records(aps,cmw,np.linspace(0,180,len(aps)).tolist())
    return [{"plan_index":1,**beam}]

def aperture_figure(cp:dict):
    left=np.array(cp["left_mm"]);right=np.array(cp["right_mm"])
    bounds=np.array(cp["leaf_boundaries_mm"])
    xjaw=cp["xjaw_mm"];yjaw=cp["yjaw_mm"]
    x0=float(min(left.min(),right.min())-12)
    x1=float(max(left.max(),right.max())+12)
    fig=go.Figure()
    fig.add_shape(type="rect",x0=x0,x1=x1,y0=bounds[0],y1=bounds[-1],fillcolor="#dff6ee",line=dict(width=0))
    for i in range(len(left)):
        fig.add_shape(type="rect",x0=x0,x1=float(left[i]),y0=float(bounds[i]),y1=float(bounds[i+1]),
                      fillcolor="#475569",line=dict(color="#cbd5e1",width=.5))
        fig.add_shape(type="rect",x0=float(right[i]),x1=x1,y0=float(bounds[i]),y1=float(bounds[i+1]),
                      fillcolor="#475569",line=dict(color="#cbd5e1",width=.5))
    if np.isfinite(xjaw).all():
        for x in xjaw: fig.add_vline(x=float(x),line_dash="dash",line_color="#dc2626")
    if np.isfinite(yjaw).all():
        for y in yjaw: fig.add_hline(y=float(y),line_dash="dash",line_color="#dc2626")
    fig.update_layout(title=f"MLC X — CP {cp['cp_index']}",xaxis_title="X, мм",yaxis_title="Y, мм",
                      height=490,margin=dict(l=15,r=15,t=50,b=25),showlegend=False)
    fig.update_yaxes(scaleanchor="x",scaleratio=1)
    return fig

st.title("RT Plan Intelligence")
st.caption("v0.1.0 · Геометрический анализ RTPLAN (MLCX) · Исследовательское ПО")
st.warning("Не использовать для клинического допуска плана. Технические показатели не заменяют дозиметрические QA-измерения.")
with st.expander("Конфиденциальность и ограничения"):
    st.markdown("""**Пациентские DICOM нельзя загружать в публичную демонстрационную версию.**
    Работайте с ними только в разрешённой защищённой локальной среде.
    Обработка выполняется в памяти процесса Streamlit; на публичном хостинге
    данные всё равно передаются на сервер. Экспорт содержит только индекс плана,
    номер луча и рассчитанные параметры, без исходного DICOM.
    Поддерживается MLCX; SAS5/SAS10 — приближённые показатели по конечным CP,
    а не временной интеграл движущейся апертуры.
    """)

# Public deployments start in synthetic-only mode: no DICOM can be uploaded.
# Explicit opt-in is required on a controlled LOCAL workstation.
local_upload = os.environ.get("RTPI_ENABLE_LOCAL_UPLOAD", "0") == "1"
a,b=st.columns([2,1])
with a:
    if local_upload:
        files=st.file_uploader("Загрузить обезличенные DICOM RTPLAN",type=["dcm","dicom"],accept_multiple_files=True)
    else:
        files=[]
        st.info("Безопасный режим: загрузка DICOM отключена. Для защищённого локального запуска установите RTPI_ENABLE_LOCAL_UPLOAD=1.")
with b: demo=st.checkbox("Синтетический пример (без DICOM)")
results=[]
if files:
    for i,file in enumerate(files,1):
        try:
            for beam in read_rtplan(BytesIO(file.getvalue()),include_control_points=True):
                results.append({"plan_index":i,**beam})
        except Exception as exc:
            st.error(f"Файл №{i}: {type(exc).__name__}: {exc}")
if demo and not files:
    results=synthetic_demo()
    st.info("Демонстрация на синтетических координатах, не на данных пациента.")

if results:
    st.subheader("Сводка по лучам/аркам")
    frame=pd.DataFrame([{k:v for k,v in r.items() if k!="control_point_data"} for r in results])
    st.dataframe(frame,hide_index=True,use_container_width=True)
    st.download_button("Скачать CSV без идентификаторов",frame.to_csv(index=False).encode("utf-8-sig"),
                       file_name="rtpi_beam_metrics.csv",mime="text/csv",type="primary")
    idx=st.selectbox("Луч/арка",list(range(len(results))),
                     format_func=lambda k:f"План №{results[k]['plan_index']} / Beam {results[k]['beam_number']}")
    result=results[idx];cps=result["control_point_data"]
    cp_index=st.slider("Контрольная точка",0,len(cps)-1,0)
    left_col,right_col=st.columns([3,2])
    with left_col: st.plotly_chart(aperture_figure(cps[cp_index]),use_container_width=True)
    with right_col:
        st.metric("Плановые MU","—" if result["monitor_units"] is None else f"{result['monitor_units']:,.1f}")
        st.metric("SAS5, MU proxy","—" if result["sas5_mu_weighted"] is None else f"{result['sas5_mu_weighted']:.3f}")
        st.metric("Leaf Travel, мм",f"{result['leaf_travel_total_mm']:,.1f}")
        st.caption("Leaf Travel учитывает все лепестки, включая вне jaws.")
    cp_data=pd.DataFrame([{k:v for k,v in cp.items() if not isinstance(v,list)} for cp in cps])
    st.subheader("Динамика контрольных точек")
    st.line_chart(cp_data.set_index("cp_index")[["aperture_area_mm2","travel_from_previous_cp_mm"]])
    st.dataframe(cp_data,hide_index=True,use_container_width=True)
else:
    st.info("Загрузите DICOM RTPLAN либо включите синтетическую демонстрацию.")

st.divider()
st.caption("MIT · Исследовательское программное обеспечение, не медицинское изделие.")
