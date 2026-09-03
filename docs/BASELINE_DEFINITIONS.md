# Baseline Definitions

Three candidate baseline definitions were evaluated because lesion-level oral dysplasia histories can be represented in more than one reasonable way.

| Definition | Field | Interpretation | Current role |
|---|---|---|---|
| Ref1 | `ref1_Lowest_Before_Thislesion` | Lowest grade before this exact lesion | Descriptive sensitivity |
| Ref2 | `ref2_Lowest_Ever_Thislesion` | Lowest grade ever for this exact lesion | Primary exploratory inference |
| Ref3 | `ref3_Lowest_AnyDx_Before_Thislesion` | Lowest diagnosis before this lesion across patient lesions | Descriptive sensitivity |

After grade filtering:

- Ref1: 22 Progressors, 1 NonProgressor
- Ref2: 25 Progressors, 3 NonProgressors
- Ref3: 25 Progressors, 2 NonProgressors

The imbalance is a central limitation and is carried through the interpretation.
