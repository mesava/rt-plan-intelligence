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
from rtpi.transitions import analyze_transitions

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
    angles=np.linspace(0,180,len(aps)).tolist()
    beam["control_point_data"]=control_point_records(aps,cmw,angles)
    trace=analyze_transitions(aps,cmw,angles,500.)
    beam["transition_summary"]=trace.summary
    beam["transition_data"]=trace.transitions
    beam["sector_data"]=trace.sectors
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
st.caption("v0.5.0 · Технические CP-переходы / VMAT-сектора · Исследовательское ПО")
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
            for beam in read_rtplan(BytesIO(file.getvalue()),include_control_points=True,
                                    include_transitions=True):
                results.append({"plan_index":i,**beam})
        except Exception as exc:
            st.error(f"Файл №{i}: {type(exc).__name__}: {exc}")
if demo and not files:
    results=synthetic_demo()
    st.info("Демонстрация на синтетических координатах, не на данных пациента.")

if results:
    st.subheader("Сводка по лучам/аркам")
    nested={"control_point_data","transition_summary","transition_data","sector_data"}
    frame=pd.DataFrame([
        {**{k:v for k,v in result.items() if k not in nested},
         "cp_max_single_leaf_step_mm":
             result.get("transition_summary",{}).get("max_single_leaf_step_mm"),
         "cp_zero_mu_leaf_motion_count":
             result.get("transition_summary",{}).get("zero_mu_leaf_motion_count")}
        for result in results
    ])
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
    trace=result.get("transition_summary",{})
    transitions=result.get("transition_data",[])
    sectors=result.get("sector_data",[])
    st.subheader("Переходы CP → CP: техническая диагностика")
    st.caption("Максимальный шаг одного лепестка, ΔMU и движение на MU. "
               "При отсутствии машинных логов это геометрия, а не фактические "
               "скорости MLC или критерий допуска плана.")
    tc=st.columns(3)
    tc[0].metric("Макс. шаг одного лепестка",
                 "—" if trace.get("max_single_leaf_step_mm") is None
                 else f"{trace['max_single_leaf_step_mm']:.2f} мм")
    tc[1].metric("Движение MLC без прироста MU",
                 str(trace.get("zero_mu_leaf_motion_count",0)))
    tc[2].metric("Поворот gantry по CP",
                 f"{trace.get('gantry_sweep_shortest_total_deg',0):.1f}°")
    if transitions:
        tdf=pd.DataFrame(transitions)
        tdf=tdf.sort_values("max_leaf_mm_per_mu",ascending=False,na_position="last")
        st.dataframe(tdf,hide_index=True,use_container_width=True)
        st.download_button("CSV по переходам CP",
                           tdf.to_csv(index=False).encode("utf-8-sig"),
                           file_name="rtpi_cp_transitions.csv",mime="text/csv")
    if sectors:
        st.subheader("Угловые сектора VMAT (по 30°)")
        st.caption("Суммы ΔMU и перемещений MLC по середине CP-интервала. "
                   "Это не вклад сектора в дозу, поэтому удаление сектора "
                   "по одной этой диаграмме недопустимо.")
        sdf=pd.DataFrame(sectors)
        st.bar_chart(sdf.set_index("sector_start_deg")["sum_delta_mu"])
        st.dataframe(sdf,hide_index=True,use_container_width=True)
        st.download_button("CSV по секторам VMAT",
                           sdf.to_csv(index=False).encode("utf-8-sig"),
                           file_name="rtpi_vmat_sectors.csv",mime="text/csv")
    else:
        st.info("Фиксированный gantry: секторный VMAT-анализ не применяется; "
                "используйте таблицу переходов CP выше.")
    st.subheader("MCSv · LSV · beamAAV (исследовательские)")
    st.caption("Jaw-aware реализация с весами ΔCMW. Эти показатели не являются подтверждением доставляемости.")
    qa_cols=st.columns(3)
    for c,key,label in zip(qa_cols,
                            ["mcs_vmat_jaw_clipped","lsv_cp_mu_weighted","aav_beam_cp_mu_weighted"],
                            ["MCSv","LSV","beamAAV"]):
        with c:
            value=result[key]
            st.metric(label,"—" if value is None else f"{value:.3f}")
    if "lsv_cp" in cp_data and "aav_beam_cp" in cp_data:
        st.line_chart(cp_data.set_index("cp_index")[["lsv_cp","aav_beam_cp"]])
    st.subheader("Динамика контрольных точек")
    st.line_chart(cp_data.set_index("cp_index")[["aperture_area_mm2","travel_from_previous_cp_mm"]])
    st.dataframe(cp_data,hide_index=True,use_container_width=True)
else:
    st.info("Загрузите DICOM RTPLAN либо включите синтетическую демонстрацию.")

st.divider()
st.caption("MIT · Исследовательское программное обеспечение, не медицинское изделие.")
