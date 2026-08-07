import gbwt.GBWT;
import gbwt.GBWTHeader;
import gbwt.bwt.BWT;
import gbwt.meta.*;
import gbwtgraph.*;
import gbwtgraph.translation.*;
import sds.RunCodec;
import sds.SdsWriter;
import java.io.*;
import java.util.*;

public class GenerateTestGBZ {
    static final int NUM_CONTIGS = 3000, NODES_PER_CONTIG = 1666;
    static final int BUBBLE_EVERY = 150, CYCLE_CONTIGS = 30, CROSS_EDGES = 100, HAPLOTYPES = 10;
    static final long SEED = 42;
    static Random rnd;
    static long nodeType(long nid) { return (nid << 1) | 0; }

    public static void main(String[] args) throws IOException {
        rnd = new Random(SEED);
        int totalNodes = NUM_CONTIGS * NODES_PER_CONTIG;
        System.out.printf("Gen: %,d contigs, %,d nodes%n", NUM_CONTIGS, totalNodes);

        List<List<long[]>> edges = new ArrayList<>(totalNodes + 1);
        edges.add(new ArrayList<>());
        for (int i=0;i<totalNodes;i++) edges.add(new ArrayList<>());

        long[] starts = new long[NUM_CONTIGS];
        String[] names = new String[NUM_CONTIGS];
        for (int c=0;c<NUM_CONTIGS;c++) { starts[c]=1L+c*NODES_PER_CONTIG; names[c]="contig_"+String.format("%04d",c); }

        boolean[] cyc = new boolean[NUM_CONTIGS];
        for (int i=0;i<CYCLE_CONTIGS;i++) cyc[rnd.nextInt(NUM_CONTIGS)]=true;

        for (int c=0;c<NUM_CONTIGS;c++) {
            long s=starts[c], e=s+NODES_PER_CONTIG-1;
            for (long n=s;n<e;n++) {
                int p=(int)(n-s);
                if (p>0 && p%BUBBLE_EVERY==0) {
                    long b1=n,b2=n+1,m=n+2; if(m>e) break;
                    edges.get((int)(n-1)).add(new long[]{nodeType(b1),0});
                    edges.get((int)(n-1)).add(new long[]{nodeType(b2),1});
                    edges.get((int)b1).add(new long[]{nodeType(m),0});
                    edges.get((int)b2).add(new long[]{nodeType(m),0});
                    n=m; continue;
                }
                if (n<e) edges.get((int)n).add(new long[]{nodeType(n+1),0});
            }
            if (cyc[c] && e>s) edges.get((int)e).add(new long[]{nodeType(s),0});
        }
        for (int i=0;i<CROSS_EDGES;i++) {
            int sc=rnd.nextInt(NUM_CONTIGS),dc=rnd.nextInt(NUM_CONTIGS);
            if (sc==dc) continue;
            long sn=starts[sc]+rnd.nextInt(NODES_PER_CONTIG), dn=starts[dc]+rnd.nextInt(NODES_PER_CONTIG);
            edges.get((int)sn).add(new long[]{nodeType(dn),0});
        }
        for (int c=0;c<NUM_CONTIGS;c++) edges.get(0).add(new long[]{nodeType(starts[c]+NODES_PER_CONTIG-1),0});
        System.out.println("  edges: "+countE(edges));

        // Build BWT records: TWO per node (forward at even index, reverse at odd)
        // Record layout: sentinel(0,1), node1_fwd(2), node1_rev(3), node2_fwd(4), node2_rev(5), ...
        int totalRecords = 2 + totalNodes * 2; // sentinel(2) + 2 per node
        List<RunCodec.NodeRecord> records = new ArrayList<>(totalRecords);
        // Sentinel: record 0 (fwd) and record 1 (rev)
        RunCodec.NodeRecord sRec = new RunCodec.NodeRecord();
        for (long[] e : edges.get(0)) sRec.outgoing.add(new long[]{e[0],0});
        int sSig = sRec.outgoing.size();
        for (int i=0;i<HAPLOTYPES;i++) sRec.runs.add(new long[]{i%sSig,1});
        records.add(sRec); // fwd
        records.add(sRec); // rev (duplicate for simplicity)

        for (int n=1; n<=totalNodes; n++) {
            RunCodec.NodeRecord fwd = new RunCodec.NodeRecord();
            List<long[]> el = edges.get(n);
            el.sort(Comparator.comparingLong(e->e[0]));
            for (long[] e : el) fwd.outgoing.add(new long[]{e[0],0});
            int sig = fwd.outgoing.size();
            if (sig>0) { int[] hc=new int[sig]; int rem=HAPLOTYPES;
                for (int e=0;e<sig-1&&rem>0;e++){hc[e]=rnd.nextInt(rem+1);rem-=hc[e];}
                hc[sig-1]=rem; for (int e=0;e<sig;e++) if(hc[e]>0) fwd.runs.add(new long[]{e,hc[e]}); }
            records.add(fwd);
            records.add(fwd); // rev = fwd for synthetic
        }
        System.out.println("  BWT records: "+records.size());

        String[] seqs = new String[totalNodes];
        char[] bases = {'A','C','G','T'};
        for (int n=0;n<totalNodes;n++) { int len=5+rnd.nextInt(5); char[] s=new char[len];
            for (int i=0;i<len;i++) s[i]=bases[rnd.nextInt(4)]; seqs[n]=new String(s); }
        System.out.println("  seqs: "+totalNodes);

        GBZFile gbz = new GBZFile();
        gbz.header = new GBZHeader(); gbz.header.tag=GBZHeader.TAG; gbz.header.version=1; gbz.header.flags=0;
        gbz.tags=new GBZTags(); gbz.tags.entries.put("source","GenerateTestGBZ");
        gbz.gbwt=new GBWT();
        gbz.gbwt.header=new GBWTHeader(); gbz.gbwt.header.tag=GBWTHeader.TAG; gbz.gbwt.header.version=5; gbz.gbwt.header.offset=1;
        gbz.gbwt.header.sequences=HAPLOTYPES*NUM_CONTIGS; gbz.gbwt.header.alphabetSize=totalNodes+1;
        gbz.gbwt.header.flags=GBWTHeader.FLAG_SIMPLE_SDS|GBWTHeader.FLAG_BIDIRECTIONAL;
        gbz.gbwt.bwt=BWT.fromNodeRecords(records);
        gbz.gbwt.header.size=gbz.gbwt.bwt.records.data.length;
        gbz.gbwt.meta=new GBWTMeta(); gbz.gbwt.meta.present=true;
        gbz.gbwt.meta.metaHeader=new MetaHeader(); gbz.gbwt.meta.metaHeader.tag=MetaHeader.TAG; gbz.gbwt.meta.metaHeader.version=2;
        gbz.gbwt.meta.metaHeader.sampleCount=1; gbz.gbwt.meta.metaHeader.haplotypeCount=HAPLOTYPES*NUM_CONTIGS;
        gbz.gbwt.meta.metaHeader.contigCount=NUM_CONTIGS; gbz.gbwt.meta.metaHeader.flags=0;
        gbz.gbwt.meta.metaSamples=new MetaSamples(); gbz.gbwt.meta.metaSamples.names=new String[]{"SAMPLE_1"};
        gbz.gbwt.meta.metaCotigs=new MetaCotigs(); gbz.gbwt.meta.metaCotigs.names=names;
        gbz.gbwt.meta.metaPaths=new MetaPaths();
        for (int c=0;c<NUM_CONTIGS;c++) for (int h=0;h<HAPLOTYPES;h++) gbz.gbwt.meta.metaPaths.paths.add(new MetaPaths.PathName(0,c,h,0));
        gbz.gbwt.header.flags|=GBWTHeader.FLAG_METADATA;
        gbz.graph=new GBWTGraph();
        gbz.graph.header=new GBWTGraphHeader(); gbz.graph.header.tag=GBWTGraphHeader.TAG; gbz.graph.header.version=3;
        gbz.graph.header.nodes=totalNodes; gbz.graph.header.flags=GBWTGraphHeader.FLAG_SIMPLE_SDS|GBWTGraphHeader.FLAG_TRANSLATION;
        gbz.graph.sequences=new GBWTGraphSequences(); gbz.graph.sequences.forward=seqs; gbz.graph.sequences.zstd=false;
        gbz.graph.translation=new Translation(); gbz.graph.translation.segments=new Segments(); gbz.graph.translation.segments.names=names;
        long[] mv=new long[NUM_CONTIGS]; for (int i=0;i<NUM_CONTIGS;i++) mv[i]=starts[i];
        gbz.graph.translation.mapping=new Mapping(); gbz.graph.translation.mapping.nodeToSegment=sds.SdsPrimitives.SparseVector.build(mv);

        String out=args.length>0?args[0]:"/Users/wenjiepeng/Downloads/test_5m.gbz";
        System.out.println("  Writing "+out+" ..."); gbz.parseTo(out);
        System.out.printf("  Done: %,d bytes%n",new File(out).length());
    }
    static long countE(List<List<long[]>> e){long c=0;for(var x:e)c+=x.size();return c;}
}
