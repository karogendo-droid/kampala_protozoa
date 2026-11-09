#!/usr/bin/env python3
"""Stage 5: purely descriptive, cross-checked reporting of archived Kampala analysis review ZIPs.

Standard library only. No FASTQ, BAM, NCBI connectivity, or fresh classifier calls.
Run from project root: python scripts/26_synthesize_stage5.py
Produces TSV tables, SVG figures, a manuscript results sketch, and a review ZIP.
"""
from __future__ import annotations
import csv, io, json, hashlib, math, sys, zipfile
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
from html import escape

ROOT=Path(__file__).resolve().parent.parent
INPUT=ROOT/'metadata'/'stage5_inputs'
OUT=ROOT/'results'/'stage5'
TAB=OUT/'tables'
FIG=OUT/'figures'
TARGETS=['Giardia','Cryptosporidium','Entamoeba','Cyclospora','Cystoisospora','Blastocystis','Enterocytozoon','Encephalitozoon']
LEADS={
 'Enterocytozoon bieneusi':'GCA_000209485.1',
 'Blastocystis ST6':'GCA_000963415.1',
 'Cyclospora cayetanensis':'GCA_002893445.1'
}
SOURCES={
 'batch':'batch_review.zip','kraken':'classification_v1_review.zip','abund':'abundance_v1_review.zip',
 'alignment':'stage3a_alignment_review.zip','blast':'stage3c_reference_specificity_review.zip',
 'audit':'stage3d_contamination_audit_review.zip','pr2':'stage3e_pr2_marker_review.zip',
 'diagnostic':'stage3f_pr2_zero_audit_review.zip','marker':'stage4_final_marker_review.zip'
}

def die(msg): raise SystemExit('STAGE 5 STOP: '+msg)

def records(source, member):
    with zipfile.ZipFile(INPUT/SOURCES[source]) as z:
        try: txt=z.read(member).decode('utf-8-sig')
        except KeyError: die(f'{SOURCES[source]} missing {member}')
    return list(csv.DictReader(io.StringIO(txt), delimiter='\t'))

def data(source, member):
    with zipfile.ZipFile(INPUT/SOURCES[source]) as z: return json.loads(z.read(member))

def integer(v, label):
    try: return int(v)
    except Exception: die(f'Non-numeric count {label}: {v!r}')

def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def comma(n):return f'{n:,}'

def svg_save(name,width,height,parts):
    txt=f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'+''.join(parts)+'</svg>'
    (FIG/name).write_text(txt,encoding='utf-8')

def svg_text(x,y,s,size=15,fill='#242424',weight='normal',anchor='start',italic=False):
    s=escape(str(s))
    return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" font-weight="{weight}" font-style="{"italic" if italic else "normal"}" text-anchor="{anchor}" fill="{fill}">{s}</text>'

def svg_rect(x,y,w,h,fill,rx=0):
    return f'<rect x="{x}" y="{y}" width="{max(0,w)}" height="{max(0,h)}" rx="{rx}" fill="{fill}"/>'

def make_fig1(samples):
    w,h=990,465
    a=[svg_rect(0,0,w,h,'#FFFFFF'),svg_text(30,38,'Figure 1. Read-pair processing by wastewater sample',21,'#222', 'bold'),
       svg_text(30,64,'Total raw, post-filter and human-screened read pairs',13,'#555')]
    colors=[('#C8C8C8','Original'),('#888888','After quality filtering'),('#B22222','After human-reference screening')]
    for j,(c,t) in enumerate(colors):
        xx=30+j*265
        a.append(svg_rect(xx,82,18,12,c));a.append(svg_text(xx+25,93,t,12))
    maxcount=max(x['raw_pairs'] for x in samples)
    plot_x,plot_w=160,785
    a.append(svg_text(plot_x,122,'0',11,'#555'))
    for k in range(0,4):
        xp=plot_x+k*plot_w/3
        a.append(f'<line x1="{xp}" y1="136" x2="{xp}" y2="418" stroke="#E8E8E8" stroke-width="1"/>')
        a.append(svg_text(xp,133,f'{maxcount*k/3/1e6:.1f}M',11,'#555',anchor='middle'))
    for i,x in enumerate(samples):
        y=155+i*43
        a.append(svg_text(30,y+17,f"{x['sample_id']} · {x['treatment_plant']}",14))
        for j,key in enumerate(['raw_pairs','qc_pairs','human_screened_pairs']):
            a.append(svg_rect(plot_x,y+j*8,plot_w*x[key]/maxcount,6,colors[j][0],rx=1))
        a.append(svg_text(952,y+18,f"{x['human_screened_pairs']/1e6:.2f}M",11,'#444',anchor='end'))
    a.append(svg_text(30,450,'Each paired-end fragment counted once. Sample-specific denominators shown in Table 1.',12,'#555'))
    svg_save('figure1_read_processing.svg',w,h,a)

def make_fig2(leads):
    w,h=900,335
    a=[svg_rect(0,0,w,h,'#FFFFFF'),svg_text(24,38,'Figure 2. Support lost after excluding flagged reference contigs',19,'#222','bold'),
       svg_text(24,62,'Pairs aligned under Stage 3A thresholds; not organism abundance',12,'#555')]
    a+=[svg_rect(28,80,14,13,'#B22222'),svg_text(48,92,'Affected by flagged contigs',12),svg_rect(270,80,14,13,'#B9B9B9'),svg_text(290,92,'Remaining, still unconfirmed',12)]
    for i,r in enumerate(leads):
        y=123+i*65
        a.append(svg_text(28,y+14,r['candidate'],15))
        bar_x=270; bar_w=500
        total=r['stage3a_pairs']; affected=r['flagged_pairs']; rem=r['remaining_pairs']
        a.append(svg_rect(bar_x,y,bar_w,25,'#EEEEEE',rx=3))
        a.append(svg_rect(bar_x,y,bar_w*affected/total,25,'#B22222',rx=3))
        a.append(svg_rect(bar_x+bar_w*affected/total,y,bar_w*rem/total,25,'#B9B9B9'))
        a.append(svg_text(bar_x+bar_w+10,y+18,f"{100*affected/total:.1f}%",13,'#222',weight='bold'))
        a.append(svg_text(270,y+43,f'{comma(affected)} flagged / {comma(rem)} remaining (total {comma(total)})',11,'#555'))
    a.append(svg_text(28,320,'Twelve selected contigs were flagged following core_nt BLAST; other reference regions were not exhaustively screened.',11,'#555'))
    svg_save('figure2_reference_flagging.svg',w,h,a)

def make_fig3(rows,controls):
    w,h=915,540
    a=[svg_rect(0,0,w,h,'#FFFFFF'),svg_text(25,39,'Figure 3. Final marker evidence across six wastewater libraries',20,'#222','bold'),
       svg_text(25,62,'Cells show isolated accepted marker read pairs; no sample–genus combination met the detection rule.',12,'#555')]
    ids=sorted({r['sample_id'] for r in rows}, key=lambda s:int(s[3:]))
    idx={(r['sample_id'],r['genus']):r for r in rows}
    x0=240;y0=127;cw=97;rh=44
    for j,s in enumerate(ids): a.append(svg_text(x0+j*cw+cw/2,103,s,13,'#242424','bold',anchor='middle'))
    for i,g in enumerate(TARGETS):
        y=y0+i*rh
        a.append(svg_text(25,y+26,g,15,italic=True))
        for j,s in enumerate(ids):
            val=int(idx[(s,g)]['accepted_pairs'])
            a.append(svg_rect(x0+j*cw+3,y+2,cw-8,rh-5,'#F3F3F3',rx=3))
            a.append(svg_text(x0+j*cw+cw/2,y+28,str(val),15,'#B22222' if val else '#999999',weight='bold' if val else 'normal',anchor='middle'))
    a.append(svg_text(25,508,'Three isolated pairs occurred in RES003 (Blastocystis), RES005 (Entamoeba), RES013 (Blastocystis).',12,'#555'))
    a.append(svg_text(25,527,'Detection requires ≥3 independent pairs and ≥2 distinct genomic positions and 100-bp reference bins.',11,'#555'))
    svg_save('figure3_marker_evidence_matrix.svg',w,h,a)

def main():
    for key,name in SOURCES.items():
        if not (INPUT/name).is_file(): die(f'Missing {INPUT/name}; extract complete Stage5 package at project root.')
    TAB.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    batch=records('batch','sample_status.tsv')
    kr=records('kraken','classification_summary.tsv')
    stage2=records('abund','genus_clade_signal_rpm.tsv')
    aln=records('alignment','all_assembly_alignment_metrics.tsv')
    blast=records('blast','reference_hotspot_blast_summary.tsv')
    flag=records('audit','affected_species_by_sample.tsv')
    refs=records('audit','flagged_reference_contigs.tsv')
    final=records('marker','stage4_genus_marker_evidence.tsv')
    ctrl=records('marker','reference_derived_control_checks.tsv')
    cov=records('marker','stage4_reference_target_coverage.tsv')
    pr2=records('pr2','pr2_marker_genus_signal_by_sample.tsv')
    diag=records('diagnostic','pr2_threshold_diagnostic_by_sample_genus.tsv')
    methods=data('marker','stage4_methods_and_limitations.json')
    ids=[r['run_accession'] for r in batch]
    if len(ids)!=6 or len(set(ids))!=6: die('Expected six unique SRA samples in batch review')
    if len(kr)!=6 or set(ids)!={r['run_accession'] for r in kr}:die('Stage1 vs batch SRA mismatch')
    if len(final)!=48 or len({(r['run_accession'],r['genus']) for r in final})!=48:die('Stage4 incomplete 6x8 matrix')
    if len({r['reference_contig'] for r in refs})!=12:die('Expected exactly 12 flagged reference contigs')
    if len(blast)!=12: die('Expected 12 BLAST hotspot summaries')
    smap={r['run_accession']:r for r in kr};samples=[]
    for x in batch:
        k=smap[x['run_accession']]
        if (x['sample_id'],x['treatment_plant'])!=(k['sample_id'],k['treatment_plant']):die('Sample metadata conflict '+x['run_accession'])
        raw=integer(x['raw_pairs'],'raw');qc=integer(x['qc_pairs'],'qc');screen=integer(x['human_screened_pairs'],'screen')
        c=integer(k['classified_pairs'],'classified');u=integer(k['unclassified_pairs'],'unclassified')
        if not raw>=qc>=screen or c+u!=screen:die('Invalid pair accounting '+x['run_accession'])
        samples.append(dict(run_accession=x['run_accession'], sample_id=x['sample_id'], treatment_plant=x['treatment_plant'],raw_pairs=raw,qc_pairs=qc,human_screened_pairs=screen,classified_pairs=c,unclassified_pairs=u,class_percent=round(100*c/screen,3),final_accepted_marker_pairs=0,final_corrob_genus_count=0))
    meta={x['run_accession']:x for x in samples}
    for r in final:
        if r['genus'] not in TARGETS or r['run_accession'] not in meta:die('Stage4 unexpected run or target')
        samp=meta[r['run_accession']]
        if (r['sample_id'],r['plant'])!=(samp['sample_id'],samp['treatment_plant']):die('Stage4 sample mismatch')
        if int(r['screened_read_pairs'])!=samp['human_screened_pairs']:die('Stage4 denominator conflict')
        n=int(r['accepted_distinct_read_pairs']); samp['final_accepted_marker_pairs']+=n
        if r['interpretation']!='NOT_CORROBORATED_AT_PRESET_THRESHOLD':samp['final_corrob_genus_count']+=1
    if sum(x['final_accepted_marker_pairs'] for x in samples)!=3:die('Unexpected isolated marker count')
    if sum(x['final_corrob_genus_count'] for x in samples):die('Expected zero corroborated genera; review changed')
    write_tsv(TAB/'table1_sample_processing.tsv',list(samples[0]),samples)
    leads=[]
    for name,acc in LEADS.items():
        aa=[r for r in aln if r['assembly_accession']==acc]
        dd=[r for r in flag if r['assembly_accession']==acc]
        if len(aa)!=6 or len(dd)!=6:die('Stage3 leading reference missing '+acc)
        pairs=sum(int(r['high_quality_distinct_pairs']) for r in aa)
        flagged=sum(int(r['pairs_with_alignment_on_flagged_contig']) for r in dd)
        remaining=sum(int(r['remaining_pairs_not_touching_flagged_contigs']) for r in dd)
        if pairs!=flagged+remaining:die('Stage3D pair reconciliation failed '+acc)
        leads.append(dict(candidate=name,assembly_accession=acc,stage3a_pairs=pairs,flagged_pairs=flagged,remaining_pairs=remaining,pct_flagged=round(100*flagged/pairs,1),stage4_genus_corroborated='No'))
    write_tsv(TAB/'table2_reference_specificity_impact.tsv',list(leads[0]),leads)
    if (sum(x['stage3a_pairs'] for x in leads),sum(x['flagged_pairs'] for x in leads))!=(9855,6251): die('Historic audit totals do not match')
    final_rows=[]
    for r in final:
        final_rows.append(dict(sample_id=r['sample_id'],plant=r['plant'],genus=r['genus'],accepted_pairs=int(r['accepted_distinct_read_pairs']),reference_positions=int(r['distinct_reference_start_sites']),reference_100bp_bins=int(r['distinct_100bp_reference_start_bins']),corroborated='No',run_accession=r['run_accession']))
    write_tsv(TAB/'table3_marker_evidence_matrix.tsv',list(final_rows[0]),final_rows)
    cover={r['genus']:int(r['pr2_ssu_reference_sequences']) for r in cov}
    genusrows=[]
    for g in TARGETS:
        f=[r for r in final_rows if r['genus']==g]
        stage2_gen=[r for r in stage2 if r['genus']==g]
        genusrows.append(dict(genus=g,pr2_reference_sequences=cover.get(g,0),exploratory_kraken_genus_clade_pairs=sum(int(r['genus_clade_pairs']) for r in stage2_gen),marker_pairs=sum(r['accepted_pairs'] for r in f),samples_with_any_marker_pair=sum(r['accepted_pairs']>0 for r in f),samples_meeting_rule=0,bugolobi_samples_meeting_rule=0,lubigi_samples_meeting_rule=0))
    write_tsv(TAB/'table4_genus_summary.tsv',list(genusrows[0]),genusrows)
    write_tsv(TAB/'table5_reference_derived_controls.tsv',list(ctrl[0]),ctrl)
    write_tsv(TAB/'supplement_pr2_diagnostic.tsv',list(diag[0]),diag)
    write_tsv(TAB/'supplement_blast_locus_findings.tsv',list(blast[0]),blast)
    make_fig1(samples);make_fig2(leads);make_fig3(final_rows,ctrl)
    raw=sum(x['raw_pairs'] for x in samples);qc=sum(x['qc_pairs'] for x in samples);screen=sum(x['human_screened_pairs'] for x in samples);assigned=sum(x['classified_pairs'] for x in samples)
    by_plant={p:[x for x in samples if x['treatment_plant']==p] for p in ['Bugolobi','Lubigi']}
    positive=[r for r in final_rows if r['accepted_pairs']]
    control_status=defaultdict(int)
    for x in ctrl:control_status[x['decision']]+=1
    report=[
      '# Stage 5. Verified synthesis of the Kampala wastewater protozoan secondary analysis','',
      'Evidence status: all target-organism assignments are exploratory; no target genus met the Stage 4 prespecified corroboration rule.','',
      '## Sequencing and classification',
      f'Six samples (three Bugolobi; three Lubigi) contained {comma(raw)} input paired-end fragments. After quality filtering {comma(qc)} pairs remained; {comma(screen)} remained after human-reference screening. Narrow-panel Kraken2 assigned {comma(assigned)} pairs ({100*assigned/screen:.2f}% of screened pairs) to any database reference. This is NOT a wastewater microbiome classification percentage.','',
      '## Reference specificity of three leading candidates',
      *[f"- {x['candidate']}: {comma(x['stage3a_pairs'])} candidate high-quality aligned pairs, {comma(x['flagged_pairs'])} associated with 12 flagged reference contigs ({x['pct_flagged']:.1f}%), {comma(x['remaining_pairs'])} remaining unconfirmed pairs." for x in leads],
      f'Together, {comma(sum(x["flagged_pairs"] for x in leads))} of {comma(sum(x["stage3a_pairs"] for x in leads))} high-quality candidate pairs ({100*sum(x["flagged_pairs"] for x in leads)/sum(x["stage3a_pairs"] for x in leads):.1f}%) were associated with flagged contigs. This is an exclusion audit of selected contigs, not a genome-wide contamination rate.','',
      '## Integrated independent marker screen',
      f'All eight target genera were evaluated in {len(final_rows)} sample–genus combinations. Zero combinations met the rule (≥3 independent pairs, ≥2 start positions and ≥2 distinct 100-bp bins). Three isolated accepted marker pairs occurred: '+', '.join(f"{x['sample_id']} {x['genus']}" for x in positive)+'.',
      f'Reference-derived control decisions: '+', '.join(f'{n} {k}' for k,n in sorted(control_status.items()))+'. These are in-silico checks, not wet-lab controls.','',
      '## Plant-level descriptive summary',
      *[f"- {plant}: {len(v)} specimens; {comma(sum(x['human_screened_pairs'] for x in v))} screened pairs; {sum(x['final_accepted_marker_pairs'] for x in v)} isolated marker pairs; zero samples meeting the genus corroboration rule." for plant,v in by_plant.items()],
      'No ecological or treatment-plant differences in protozoan occurrence can be inferred.','',
      '## Reporting boundaries',
      '- Do not claim parasite absence, species/subtype confirmation, organism abundance, viability, or host infection.',
      '- Only 12 selected reference contigs were flagged; contamination of all reference sequences was not exhaustively assessed.',
      '- The narrow bacterial comparator panel makes the Kraken2 classified fraction nonrepresentative of the microbial community.',
      '- The marker method lacks empirically established limits of detection and extraction/PCR controls; ambiguous in-silico controls for Cyclospora and Cystoisospora limit interpretation.',
      '- The original proposal included bacterial-indicator concordance, which has NOT been performed; do not claim this objective was met.',
      '- Sampling source/mode and time differ; three specimens per plant, descriptive comparison only. No seasonality/treatment-removal or population prevalence claims.',
      '- User should verify original specimen metadata and ethics requirements before manuscript submission.','',
      '## Reproducibility',
      'Input ZIP hashes are recorded in stage5_manifest.json. The Stage5 script operates only on archived aggregate reports. No FASTQ/BAM read extraction is performed.',''
    ]
    (OUT/'stage5_results_summary.md').write_text('\n'.join(report),encoding='utf-8')
    mani=dict(stage='5',generated_utc=datetime.now(timezone.utc).isoformat(),source_sha256={name:hashlib.sha256((INPUT/name).read_bytes()).hexdigest() for name in SOURCES.values()},sample_count=len(samples),target_genera=TARGETS,screened_pairs=screen,classified_pairs=assigned,marker_supported_pairs=sum(x['accepted_pairs'] for x in final_rows),sample_genus_combinations=len(final_rows),corroborated_combinations=0,controls=dict(control_status),stages=['batch','1','2','3A','3C','3D','3E','3F','4'],limitations=methods.get('caveats',[]))
    (OUT/'stage5_manifest.json').write_text(json.dumps(mani,indent=2)+'\n',encoding='utf-8')
    package=ROOT/'metadata'/'taxonomy'/'stage5_results_review.zip';package.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as zz:
        for p in [*TAB.glob('*.tsv'), *FIG.glob('*.svg'),OUT/'stage5_results_summary.md',OUT/'stage5_manifest.json']:
            zz.write(p,p.relative_to(ROOT))
    print('STAGE 5 COMPLETE: all 9 source reviews reconciled')
    print(f'Six samples: {comma(screen)} screened pairs; {comma(assigned)} narrow-panel Kraken2 classified')
    print(f'Flagged leading-candidate pairs: {comma(sum(x["flagged_pairs"] for x in leads))}/{comma(sum(x["stage3a_pairs"] for x in leads))}')
    print(f'Isolated marker pairs: {len(positive)}; corroborated target sample–genus combinations: 0/48')
    print('Output tables / figures:',OUT)
    print('Review ZIP:',package)

if __name__=='__main__':main()
