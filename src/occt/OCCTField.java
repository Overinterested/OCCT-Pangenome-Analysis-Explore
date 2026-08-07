package occt;

import edu.sysu.pmglab.ecc.field.FieldSchema;
import edu.sysu.pmglab.ecc.field.IFieldSchema;
import edu.sysu.pmglab.ecc.type.FieldType;

/**
 * @author Wenjie Peng
 * @create 2026-08-04 16:52
 * @description
 */
public enum OCCTField {
    // 一个示例，需要把OCCTField分为多个group，每个group含有多个子字段
    A_GROUP(
            "X",
            new String[]{"X_A", "X_B", "X_C"},
            new FieldType[]{FieldType.int32List, FieldType.int32, FieldType.float32}
    ),
    B_GROUP(
            "XX",
            new String[]{"XX_A", "XX_B", "XX_C"},
            new FieldType[]{FieldType.bytecode, FieldType.string, FieldType.stringList}
    ),;

    final String group;
    final IFieldSchema metas;
    final String[] groupNames;
    final FieldType[] groupFields;

    OCCTField(String group, String[] groupNames, FieldType[] groupFields) {
        this.group = group;
        this.groupNames = groupNames;
        this.groupFields = groupFields;
        FieldSchema fieldMetas = new FieldSchema();
        for (int i = 0; i < groupNames.length; i++) {
            fieldMetas.addField(group, groupNames[i], groupFields[i]);
        }
        this.metas = fieldMetas;
    }

    public IFieldSchema getMetas() {
        return metas;
    }

}
