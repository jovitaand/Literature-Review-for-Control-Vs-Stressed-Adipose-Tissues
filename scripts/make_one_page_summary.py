import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
D="/usr/share/fonts/truetype/dejavu/"
pdfmetrics.registerFont(TTFont("DV",D+"DejaVuSans.ttf")); pdfmetrics.registerFont(TTFont("DV-B",D+"DejaVuSans-Bold.ttf"))
pdfmetrics.registerFontFamily("DV",normal="DV",bold="DV-B",italic="DV",boldItalic="DV-B")
b=pd.read_csv("GSE268904_barcodes_clean.csv",usecols=["sample","condition","endo_specific","total_umi"])
order={"lean":0,"obese":1,"obese diabetic":2}; nice={"lean":"Lean","obese":"Obese, no diabetes","obese diabetic":"Obese with type 2 diabetes"}
col={"lean":"#2E75B6","obese":"#7F8FA6","obese diabetic":"#C55A11"}
g=b.groupby(["sample","condition"],as_index=False).agg(n=("total_umi","size"),e=("endo_specific","sum")); g["pct"]=100*g.e/g.n
g["o"]=g.condition.map(order); g=g.sort_values(["o","pct"],ascending=[True,False]).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(7.6,3.7),dpi=200)
ax.bar(range(len(g)),g.pct,color=[col[c] for c in g.condition],width=0.72)
for i,r in g.iterrows(): ax.text(i,r.pct+0.8,f"{r.pct:.0f}%",ha="center",fontsize=7)
start=0
for c in ["lean","obese","obese diabetic"]:
    idx=g.index[g.condition==c]; m=g.pct[idx].mean()
    ax.hlines(m,idx.min()-0.4,idx.max()+0.4,colors="black",linestyles="--",lw=1.1)
    ax.text((idx.min()+idx.max())/2,-8.5,f"{nice[c]}\n(average {m:.0f}%, {len(idx)} people)",ha="center",va="top",fontsize=7.5,color=col[c],fontweight="bold")
ax.set_xticks(range(len(g))); ax.set_xticklabels(g["sample"],fontsize=7,rotation=0)
ax.set_ylabel("Estimated blood-vessel cells\n(% of all cells measured)",fontsize=8); ax.set_ylim(0,48)
for s in ["top","right"]: ax.spines[s].set_visible(False)
plt.tight_layout(rect=(0,0.08,1,1)); fig.savefig("chart.png"); plt.close(fig)
tot=len(b); ev=int(b.endo_specific.sum())
S=lambda **k: ParagraphStyle("s",**{"fontName":"DV",**k})
h1=S(fontSize=17,leading=21,textColor=colors.HexColor("#1F4E79"),fontName="DV-B"); sub=S(fontSize=9,leading=12,textColor=colors.HexColor("#595959"))
h2=S(fontSize=10.5,leading=13,textColor=colors.HexColor("#1F4E79"),fontName="DV-B",spaceBefore=6,spaceAfter=2); body=S(fontSize=8.6,leading=11.6); small=S(fontSize=7.4,leading=9.6,textColor=colors.HexColor("#595959"))
bnum=S(fontSize=18,leading=20,textColor=colors.HexColor("#1F4E79"),fontName="DV-B",alignment=1); blab=S(fontSize=7.8,leading=9.5,alignment=1,textColor=colors.HexColor("#595959"))
doc=SimpleDocTemplate("GSE268904_one_page_summary.pdf",pagesize=A4,leftMargin=16*mm,rightMargin=16*mm,topMargin=12*mm,bottomMargin=10*mm)
st=[Paragraph("Blood-vessel cells in human fat: lean, obese, and obese with type 2 diabetes",h1),
 Paragraph("Public dataset GSE268904: fat tissue under the skin from 14 people. Counts below are estimates from gene activity, not a laboratory identification.",sub),Spacer(1,5)]
kp=Table([[Paragraph(f"{tot:,}",bnum),Paragraph("14",bnum),Paragraph(f"{ev:,}",bnum),Paragraph(f"{100*ev/tot:.0f}%",bnum)],
 [Paragraph("individual cells (nuclei) measured",blab),Paragraph("people (samples)",blab),Paragraph("estimated blood-vessel cells",blab),Paragraph("of all cells",blab)]],colWidths=[44*mm]*4)
kp.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#F2F6FA")),("VALIGN",(0,0),(-1,-1),"MIDDLE"),("BOTTOMPADDING",(0,1),(-1,1),5),("TOPPADDING",(0,0),(-1,0),5)]))
st+= [kp,Spacer(1,4),Image("chart.png",width=176*mm,height=176*mm*3.7/7.6)]
h=b.groupby("condition").agg(p=("sample","nunique"),n=("total_umi","size"),e=("endo_specific","sum")); gp=g.groupby("condition").pct.agg(["mean","min","max"])
rows=[["Group","People","Cells measured","Est. blood-vessel cells","Average per person","Range across people"]]
for c in ["lean","obese","obese diabetic"]: rows.append([nice[c],int(h.p[c]),f"{int(h.n[c]):,}",f"{int(h.e[c]):,}",f"{gp['mean'][c]:.0f}%",f"{gp['min'][c]:.0f}% to {gp['max'][c]:.0f}%"])
t=Table([[Paragraph(str(x),small if i else S(fontSize=7.4,leading=9.6,textColor=colors.white,fontName="DV-B")) for x in r] for i,r in enumerate(rows)],colWidths=[44*mm,17*mm,27*mm,32*mm,29*mm,27*mm])
t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#1F4E79")),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#F2F6FA")]),("GRID",(0,0),(-1,-1),0.25,colors.HexColor("#D0DCE8")),("TOPPADDING",(0,0),(-1,-1),2.5),("BOTTOMPADDING",(0,0),(-1,-1),2.5)]))
st+=[Spacer(1,3),t,
 Paragraph("What this shows",h2),
 Paragraph("About 1 in 5 cells measured looks like a blood-vessel lining cell. The share is very different from person to person (from "+f"{g.pct.min():.0f}% to {g.pct.max():.0f}%"+"), and the groups overlap heavily.",body),
 Paragraph("Please read with care",h2),
 Paragraph("&#8226; Only 3 people have type 2 diabetes, so group differences are exploratory and prove nothing yet.<br/>&#8226; The obese groups are much heavier (BMI about 50) than the lean group (about 23), so any difference may reflect weight, not diabetes.<br/>&#8226; Cells are not people: the real sample size is 14.",body),
 Paragraph("How a cell was counted as a blood-vessel cell (<b>endo_specific</b>)",h2),
 Paragraph("A cell is counted when (1) at least 3 of 6 blood-vessel marker genes (PECAM1, VWF, CLDN5, EGFL7, FLT1, PTPRB) are active in it, (2) blood-vessel genes are its strongest signal compared with fat, connective-tissue, immune and vessel-wall genes, and (3) that signal is at least twice as strong as any of the others. This guards against cells that only look blood-vessel-like because of stray RNA in the sample. It is still a marker-based estimate, not a formal cell-type annotation.",body),
 Spacer(1,3),Paragraph("Source: GEO GSE268904 (Aarhus University), analysed with scripts/GSE268904_extract.R. Counts: barcodes_clean table, 131,641 nuclei.",small)]
doc.build(st)
print("built")
