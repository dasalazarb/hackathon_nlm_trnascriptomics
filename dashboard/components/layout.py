import streamlit as st


def intro(question: str, text: str) -> None:
    st.markdown(f"### {question}")
    st.write(text)


def source_box(*names: str) -> None:
    with st.expander("Evidence / Source"):
        st.write("R pipeline outputs used in this view:")
        for name in names:
            st.code(name)
