import json, re, sys
from pathlib import Path
SP = Path('/tmp/claude-1000/-home-chen-project-Immrep26/1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/deck')
OUT = Path('/home/chen/project_Immrep26/reports/IMMREP26_讲解.html')
deck = json.loads((SP/'project/deck.json').read_text())
FB = {"'Noto Sans SC', Arial, sans-serif":"'Noto Sans SC','PingFang SC','Hiragino Sans GB','Microsoft YaHei',sans-serif",
      "'Noto Serif SC', Georgia, serif":"'Noto Serif SC','Songti SC','SimSun',Georgia,serif",
      "'JetBrains Mono', 'Courier New', monospace":"'JetBrains Mono','SFMono-Regular',Consolas,'Courier New',monospace"}
secs, notes = [], []
for sid in deck['order']:
    h = (SP/f'project/slides/{sid}.html').read_text()
    for a,b in FB.items(): h = h.replace(a,b)
    m = re.search(r'<aside>(.*?)</aside>', h, re.S)
    notes.append((sid, m.group(1).strip() if m else ''))
    secs.append(re.sub(r'<aside>.*?</aside>','',h,flags=re.S).strip())
css = """*{box-sizing:border-box}
html,body{margin:0;padding:0;background:#2a3133}
/* The slide runtime treats margin as a no-op; a browser does not, and the
   default margins on p/h2 add ~60px per block. */
p,h1,h2,h3,ul,ol,table{margin:0}
section{width:1920px;height:1080px;overflow:hidden;position:relative;
  display:flex;flex-direction:column;margin:0 auto 24px;
  box-shadow:0 2px 24px rgba(0,0,0,.35)}
table{border-collapse:collapse;width:100%}
th,td{padding:7px 14px 7px 0;text-align:left;vertical-align:baseline}
sup{font-size:.65em;vertical-align:super}
/* <x-shape> is the slide runtime's own element; a plain browser renders it as
   a filled box. Clip it into the arrow the runtime would have drawn, so the
   standalone file and the runtime show the same diagram. */
x-shape{display:inline-block;vertical-align:middle;flex:0 0 auto}
x-shape[kind="arrow-right"]{clip-path:polygon(0 32%,58% 32%,58% 0,100% 50%,58% 100%,58% 68%,0 68%)}
x-shape[kind="arrow-down"]{clip-path:polygon(32% 0,32% 58%,0 58%,50% 100%,100% 58%,68% 58%,68% 0)}
@media screen{body{padding:24px 0}.stage{display:flex;flex-direction:column;align-items:center}}
@page{size:1920px 1080px;margin:0}
@media print{html,body{background:#fff}body{padding:0}.stage{zoom:1 !important}
  section{margin:0;box-shadow:none;break-after:page;page-break-after:always}
  section:last-of-type{break-after:auto;page-break-after:auto}
  .notes{display:none !important}}"""
js = """(function(){var st=document.querySelector('.stage');
function fit(){st.style.zoom=Math.min(1,(window.innerWidth-36)/1920);}
if(!window.matchMedia('print').matches){fit();window.addEventListener('resize',fit);}})();"""
faces = "\n".join(f'<link rel="stylesheet" href="{f["href"]}">' for f in deck['faces'].values())
nh = "\n".join(
  f'<div style="max-width:1080px;margin:0 auto 20px;padding:18px 24px;background:#fbfcfc;'
  f'border-left:3px solid #15635c;font-family:\'Noto Sans SC\',\'PingFang SC\',sans-serif;'
  f'font-size:15px;line-height:1.8;color:#11181b">'
  f'<b style="font-family:monospace;color:#15635c">{i+1:02d} · {sid}</b><br>{n}</div>'
  for i,(sid,n) in enumerate(notes) if n)
doc = f"""<!doctype html>
<html lang="zh-Hans"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{deck['title']}</title>
{faces}
<style>{css}</style></head><body>
<div class="stage">
{chr(10).join(secs)}
</div>
<div class="notes" style="background:#2a3133;padding:36px 20px 56px">
<h2 style="max-width:1080px;margin:0 auto 20px;font-family:'Noto Serif SC',serif;color:#e7edeb;font-size:24px">演讲备注</h2>
{nh}</div>
<script>{js}</script></body></html>"""
OUT.write_text(doc, encoding='utf-8')
probe = """<script>window.addEventListener('load',function(){
document.querySelector('.stage').style.zoom='';var out=[];
document.querySelectorAll('section').forEach(function(el){var inner=0;
el.childNodes.forEach(function(c){if(c.getBoundingClientRect){var r=c.getBoundingClientRect();inner=Math.max(inner,r.bottom);}});
var limit=1080-parseFloat(getComputedStyle(el).paddingBottom);
var used=Math.round(inner-el.getBoundingClientRect().top);
if(used>limit)out.push(el.id+' '+used+'>'+Math.round(limit));});
document.title='OVER '+(out.join(' | ')||'NONE');});</script>"""
Path('/mnt/c/Users/Public/chk.html').write_text(doc.replace('</body>', probe+'</body>'), encoding='utf-8')
print(f"{len(secs)} slides, {OUT.stat().st_size:,} bytes")
