def test_load_module():
    from site_annotate.utils import data_generator
    data_generator.generate_test_data # this will fail if it doesn't exist. so if it passes it exists


from site_annotate.utils import data_generator


def test_generate_data():
    df = data_generator.generate_test_data(1)
    assert "Peptide Length" in df.columns
