# tests/test_stratify_shortcodes.py
import pandas as pd
from CRISPRSCope.cli import stratify_data

def test_stratify_emits_only_short_codes():
    df = pd.DataFrame({
        'Barcode': ['BC1','BC2','BC3','BC4'],
        'Barcode Rank': [1, 500, 15000, 20000],
        'Amplicon Score': [2.0, 0.2, 1.5, 0.0],
    }).set_index('Barcode')

    out = stratify_data(df.copy())
    # must have a Color column
    assert 'Color' in out.columns

    # only allowed short codes appear
    allowed = {'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}
    seen = set(out['Color'].unique().tolist())
    bad = seen - allowed
    assert not bad, f"Unexpected codes present: {bad}"

    # spot-check expected assignments from canonical thresholds
    # BC1: rank=1, score=2.0 => HQ_HI
    assert out.loc['BC1','Color'] == 'HQ_HI'
    # BC2: rank=500, score=0.2 => LQ_HI
    assert out.loc['BC2','Color'] == 'LQ_HI'
    # BC3: rank=15000, score=1.5 => HQ_LO
    assert out.loc['BC3','Color'] == 'HQ_LO'
    # BC4: rank=20000, score=0.0 => LQ_LO
    assert out.loc['BC4','Color'] == 'LQ_LO'
