/**
 * {@link GBZFile#parse(String, GBZFilter)} 的模块过滤开关：
 * 置为 false 的模块在读取时直接跳过（只读长度字段，数据部分 seek 过去），
 * 不在内存里实例化。被跳过的模块在返回的 GBZFile 里保持为空壳对象
 * （例如 meta.present == false、bwt.index == null），请不要使用其内容。
 *
 * 注意：过滤后的 GBZFile 用于分析，不适合再 parseTo() 写回——
 * 缺了 BWT/图部分的文件不是有效的 GBZ（仅跳过 metadata/DASamples 时写回仍有效，
 * 因为 GBWT.encode 会同步清掉 header 里对应的 flag）。
 */
public class GBZFilter {

    public boolean gbzTags = true;
    public boolean gbwtTags = true;
    public boolean bwt = true;
    public boolean daSamples = true;
    public boolean metadata = true;
    public boolean graphSequences = true;
    public boolean graphTranslation = true;

    public static final GBZFilter LOAD_ALL =  new GBZFilter();

    /** 全部模块都加载。 */
    public static GBZFilter newInstance() { return new GBZFilter(); }

    /** 只要 header/tags/metadata，跳过 BWT、DASamples 和整个图部分（文件中最大的几块）。 */
    public static GBZFilter metadataOnly() {
        return newInstance().withoutBWT().withoutDASamples()
                .withoutGraphSequences().withoutGraphTranslation();
    }

    /** 示例场景：不保存 GBWT 的 metadata 部分。 */
    public GBZFilter withoutMetadata() { this.metadata = false; return this; }
    public GBZFilter withoutDASamples() { this.daSamples = false; return this; }
    public GBZFilter withoutBWT() { this.bwt = false; return this; }
    public GBZFilter withoutGBZTags() { this.gbzTags = false; return this; }
    public GBZFilter withoutGBWTTags() { this.gbwtTags = false; return this; }
    public GBZFilter withoutGraphSequences() { this.graphSequences = false; return this; }
    public GBZFilter withoutGraphTranslation() { this.graphTranslation = false; return this; }
}
