import requests, time, json, re
DOIS={"Bhardwaj2002":"10.1124/jpet.102.034728","Buckley2018":"10.1126/scitranslmed.aar7047","Collnot2006":"10.1016/j.jconrel.2005.11.005",
"Dintaman1999":"10.1023/a:1015000503629","Griesser2017":"10.1016/j.ijpharm.2017.02.019","Maher2016":"10.1016/j.addr.2016.06.005","Ristroph2019":"10.1039/c9na00308h",
"Twarog2019":"10.3390/pharmaceutics11020078","Overgaard2021":"10.1007/s40262-021-01025-x","Granhall2019":"10.1007/s40262-018-0728-4","Brayden2020":"10.1016/j.addr.2020.05.007",
"Pechenov2021":"10.1038/s41598-021-01750-0","Drucker2020":"10.1038/s41573-019-0053-0","Atuma2001":"10.1152/ajpgi.2001.280.5.g922","Cu2009":"10.1016/j.addr.2008.09.006",
"Saltelli2010":"10.1016/j.cpc.2009.09.018","Sobol2001":"10.1016/s0378-4754(00)00270-6","Jansen1999":"10.1016/s0010-4655(98)00154-4","Iman1982":"10.1080/03610918208812265",
"Beaumont2002":"10.1093/genetics/162.4.2025","Marino2008":"10.1016/j.jtbi.2008.04.011"}
out={}
for k,d in DOIS.items():
    for t in range(4):
        try:
            m=requests.get("https://api.crossref.org/works/"+d,params={"mailto":"sandler.leon@gmail.com"},timeout=60).json()['message']; break
        except Exception: time.sleep(3)
    au=m.get('author',[]); 
    def ini(g): return "".join(x[0]+"." for x in re.split(r"[\s-]+",g) if x) if g else ""
    names=[f"{a.get('family','')}, {ini(a.get('given',''))}" for a in au]
    year=(m.get('issued',{}).get('date-parts') or [[None]])[0][0]
    title=re.sub('<[^>]+>','',m['title'][0]).strip().rstrip('.')
    cont=(m.get('container-title') or [''])[0]; vol=m.get('volume',''); iss=m.get('issue',''); pg=m.get('page','') or m.get('article-number','')
    out[k]=dict(authors=names,year=year,title=title,journal=cont,volume=vol,issue=iss,pages=pg,doi=d)
    print(k,len(names),year,cont[:30],vol,iss,pg)
    time.sleep(1)
json.dump(out,open('build/refs_meta.json','w',encoding='utf-8'),indent=1,ensure_ascii=False)
