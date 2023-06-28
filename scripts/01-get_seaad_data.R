
#get all SEAAD data
source('../../utils/r_utils.R')
library(Seurat)
out<-'../singlecell-brain/outputs/01-SEAAD_data'
dir.create(out,recursive = T)

# system('curl -o ref-data/SEAAD/oligo.rds "https://corpora-data-prod.s3.amazonaws.com/0bcc5dbd-0a6b-461a-b10c-d79a326a29cf/local.rds?AWSAccessKeyId=ASIATLYQ5N5XTJULZ6FX&Signature=l1L7VKdhl9PUA%2BoGc0N6FAW6Q6w%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIHANI3tq2E2kkqgt9rN9bFItFqqLjG%2F%2BukVNlqm21qr0AiA6wTlG6NhNqvyVqxQ9AEETFBhyKfNvzxei0yWNJ7ah1Sr0Awj0%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAEaDDIzMTQyNjg0NjU3NSIMJlE6hpbBIxFtFYkyKsgDtquTcniXhHXdRKg%2FD6oVy15CNAm9ad16qfpdPrmuv%2B8cnhVZGA0%2FVeWvza49oG4Wgj%2B%2BBt5T7wsKWdJbPbrfIDIkb2C5c93B0LsRAxDemoGOtJ7uyD3Gz%2FIqyD7pspOxGtuQdIaeMMtaTOXsumn3Q9wR%2FrmofTnbbx3Re3oeKmemb27QYTGHJxFKflx%2B%2BJuQLn6OL%2Bz1GeNRtaCq9g2iLfkV7lWzqHBM5MiqpuMRIF1bPSbjDFVAkK8szfmHVC9JO%2BtBJsOiAjauXTYIGLYWBBYl%2B%2FOhjrvwIxOkyzQDcFPXjO2Xi0QBdSmSWjWswTz8dZEKggap2vOqea4mnRPNU3ZLAHhWPoZfokcvY4MImqpVWza%2BJO5BRsutgUCiwgMWqE%2FM5PhIigACj1cZIo%2F1pRvhTK7wcgCxotRpZ1EqPn4ctgLH4tNUbpwulI4S%2Fbpmr%2FpHImpJ9hLYgQv5XiSiKAumGelOnqNdW3qVFk5XxR%2Ft6aSTCxQs2zmbCyZOk9Kd6lOai4LgaR1tHUmrH%2FB1DlrMh%2B7g3Bpe7bRHHippTRx7c%2BXqSxF7317XOhXkyeHWIObnOtmFBoiqlybKNq72RwuGwIcsdBCsMMbGuaMGOqYBo0dUCzf6veK4K7I%2FZPdhwR1nglVfOnGxCOwQd3na4ks%2F4gFRPd%2FLoEndEgSaip%2FSA%2FiNfIngGL34yiadYo9eJ9xXxxK%2FUJ0eP%2BCAdn44PTCw893lg5fXWhVU5BoK0cwbzybnZKcctlrL4XxVVcrgiy0Kv8NO5y6Ck9Qy4S0LRLvhUm0k45qSnwqQdOpMQsv0u7yHmXxhwBn9AfIOcom8eTO4I240Lg%3D%3D&Expires=1685567849"')
# system('curl -o ref-data/SEAAD/Vip.rds "https://corpora-data-prod.s3.amazonaws.com/14632982-a230-40f3-9873-600d4ddec5aa/local.rds?AWSAccessKeyId=ASIATLYQ5N5XV2SDXEL5&Signature=FZwItgv7pcbheWzY9lpHxF1xIP8%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCID3h3jAtF8yhD5RZrNPzmbhbWTfErrDZfhoSWWkeIwCcAiEAzAMU0FhETPpyBAebZTV%2BjbCcGav7T7%2BeYO%2B38iZIbMkq9AMI8%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgwyMzE0MjY4NDY1NzUiDPATOnyPkr085qmBzCrIA3Bwhue0%2FaxyfGFqOE1Y2IFqZiOYQ65p73TFMHsr%2FdnCQwehzEA%2BsWxCpieiyWCZBN61%2FjxE%2FT3xhVGIj16QyGQqiLAkuWPkvNmp9o60urX5po59WgU1CyX9Gp8Kx6rra3rhE5ckzmUlFU4vkiClPZEc6baQKVxqPIuLYynFlnaDXwtWZdvz35nGVkGVLSG53yyTFhb79GSFaMuTDosZS%2BZ686dtRUaxGCTbRm%2FNBhIU1LEmp2WkB9qpo5Tj3lRC7OEXSN0%2FwGEWhCMgrhNdLFq9Xi342ZqjO7SkcSSQE85u5EmY2kgl887tTzmV4QCQx4Bm6AcSZBes%2FneWWSBboaZRHgs%2Fnp29n%2F8BCzp%2BMpeE8pgwHNwV%2BzS%2Bm5d3u72%2FOBqNjtEveuLU5cGhLnxBjs5rNOAK%2Bu0XwC%2BsJbz3Ah42vsVqPjW%2FoXc6fUBZsOF6C15BCXmIqSKeTQHLXCdNTHgC0nJfumMUJio99GBhn9uIGB5x04hcX9gp6YRoUk2IEI2Oj3XaGC7nWZLNtPxUGSBJ5ne8suEtpFlaRZMIRUWbhns2uf7T4rbLVWIb1j%2BJdbsww6FgSUmLJ8aB0WPTqSR2SSsjj3BTkzCeqrmjBjqlAevXbGwDpQ03KR6r02VfrEuPhw%2B98UqZuKFh9%2BX13Ootdr4XGWMkK0HRkqMlij2zsla%2BkQL4TZjgANz3n6Q4nAAry7dKfDtZsm06Luc6U1sEJWcbMG19Q08BTocCFr5uEKkXqKPYklCoApXjldV3Sm6odP8A4EvAsUbB9%2BP6bOlpkDlPwCBrVIiqDdVRgVD9bCGPErJ77GYAy2NJy84lVuOyw2q%2B7w%3D%3D&Expires=1685569296"')
# system('curl -o ref-data/SEAAD/Pvalb.rds "https://corpora-data-prod.s3.amazonaws.com/705feb0c-0733-41a7-8887-ad88db62464a/local.rds?AWSAccessKeyId=ASIATLYQ5N5XTJULZ6FX&Signature=T9PpIfeAb6lo1H8NPlP30qpby7I%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIHANI3tq2E2kkqgt9rN9bFItFqqLjG%2F%2BukVNlqm21qr0AiA6wTlG6NhNqvyVqxQ9AEETFBhyKfNvzxei0yWNJ7ah1Sr0Awj0%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAEaDDIzMTQyNjg0NjU3NSIMJlE6hpbBIxFtFYkyKsgDtquTcniXhHXdRKg%2FD6oVy15CNAm9ad16qfpdPrmuv%2B8cnhVZGA0%2FVeWvza49oG4Wgj%2B%2BBt5T7wsKWdJbPbrfIDIkb2C5c93B0LsRAxDemoGOtJ7uyD3Gz%2FIqyD7pspOxGtuQdIaeMMtaTOXsumn3Q9wR%2FrmofTnbbx3Re3oeKmemb27QYTGHJxFKflx%2B%2BJuQLn6OL%2Bz1GeNRtaCq9g2iLfkV7lWzqHBM5MiqpuMRIF1bPSbjDFVAkK8szfmHVC9JO%2BtBJsOiAjauXTYIGLYWBBYl%2B%2FOhjrvwIxOkyzQDcFPXjO2Xi0QBdSmSWjWswTz8dZEKggap2vOqea4mnRPNU3ZLAHhWPoZfokcvY4MImqpVWza%2BJO5BRsutgUCiwgMWqE%2FM5PhIigACj1cZIo%2F1pRvhTK7wcgCxotRpZ1EqPn4ctgLH4tNUbpwulI4S%2Fbpmr%2FpHImpJ9hLYgQv5XiSiKAumGelOnqNdW3qVFk5XxR%2Ft6aSTCxQs2zmbCyZOk9Kd6lOai4LgaR1tHUmrH%2FB1DlrMh%2B7g3Bpe7bRHHippTRx7c%2BXqSxF7317XOhXkyeHWIObnOtmFBoiqlybKNq72RwuGwIcsdBCsMMbGuaMGOqYBo0dUCzf6veK4K7I%2FZPdhwR1nglVfOnGxCOwQd3na4ks%2F4gFRPd%2FLoEndEgSaip%2FSA%2FiNfIngGL34yiadYo9eJ9xXxxK%2FUJ0eP%2BCAdn44PTCw893lg5fXWhVU5BoK0cwbzybnZKcctlrL4XxVVcrgiy0Kv8NO5y6Ck9Qy4S0LRLvhUm0k45qSnwqQdOpMQsv0u7yHmXxhwBn9AfIOcom8eTO4I240Lg%3D%3D&Expires=1685569333"')
# system('curl -o ref-data/SEAAD/Sst.rds "https://corpora-data-prod.s3.amazonaws.com/673ff65f-a98b-49f3-a4a9-175f93a36860/local.rds?AWSAccessKeyId=ASIATLYQ5N5X7ILTHZSG&Signature=PCmPQwGhq3gFM0D9aseUX9wTqVk%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCmk%2BDKAO6TSaqlAMyp8xuGfdAUo6xyoyQX2FL0QymLcwIhAMdmF%2FPA1SQ4zAJg7nmsHlZ%2F9gjda274if1h0dexdS%2BfKvQDCPX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igz9aRgBioj4Xd7Gd70qyAOU13Cq8CtyTVdeynRyHPtB8bMBkLty%2B5dIVx325rmT5yL5D9NmaB%2BOpvuwMpID%2FmQGm3RjCb7B3wnEwu7x%2Bcu%2BhkBVM%2BPBfdBc6jhdqW8X0nvA%2BmhgJF6qA2GtiI%2BuHfapUyMTKupUrbGd%2Fr8wTNjqZj%2B%2BHLs11Mr2vIlIIh0iK2NrMr7BHHpjew25iVfDjoJmvpmnJCI7mLZuvSkXLFhM5SDCMkVscQ%2BX1BrTBV09ImP9uxQKz%2B5AIAUHZaqNxWGrguieesqWFTwsRBhuUyONegDF9nSipNrHpj0G4OlcGdO6pyM2yvcenAfxONM7pVpSl6kZ9SGFGbZU91ggjfSuu3nrfHm%2B11QEuNCn62MPru4il%2FrWAOjy2lYti6j%2B17I2j19m5fw3%2FJcsTuF8f9ZIN0wEMUMMLpWAYFHALV%2FGDVmu60WH2cAKVITeCOkdp7d7r%2BK4KQw7yg8PgPsIJlfTRoBK5rywiUuV9z8VfmioGVjqGkLAN0DL6xSIopTGvOSnlgFIjwf0mRa9OV5mHOb1DsLP5UElXmlA4%2BpRZKhQL%2BKxupPxn1SHZWS6sGEeeVYF0pls%2FRfIPO3u5wlkqv1Iq9wPg4JNd4Mw%2Bse5owY6pAFuZMwAir42WfJK0XtVwhqPRxdzcsNGdFHypq%2BOYwbbOrNRnqxg3A13kjpPaY3%2B118HmaMA2YbqLb9%2BWBDHjyD7iRr2SjX%2BMLcAa0UxnOP%2FFWAdSIv2kvuuyb4W0kuCXmAZNkIPdOwyO%2FCyU3spBg40jDpnlJTFu%2BsFmz3if3yFIO4KL%2Fed3WGxXjqBcdDMAhs5D9sNwZsNg6hIUD9GqFybPHsC3w%3D%3D&Expires=1685569368"')
# system('curl -o ref-data/SEAAD/Lamp5.rds "https://corpora-data-prod.s3.amazonaws.com/108ac351-76d3-46c0-8d82-8a45a3e244b3/local.rds?AWSAccessKeyId=ASIATLYQ5N5X572JEX7C&Signature=7iKHZ8nT1dBvRkTMf7zFKHBRO1Y%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMr%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCSLYD6c9GJDJ462R3fyKA%2BtqLAfd3ZLYdTs2l8Sm8OqAIhAOufFKbvnALqcS51UPJRams4iUvBfeoo4DEazM%2BCELe8KvQDCPP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1IgzVqIuk562TH8S9imEqyAN32f3C0LsQrgnOZHSbd%2F35aLvFO88avQ3fpI%2B9aFN4d4cQhj4%2FAXOaGvo1pTcjmSZvzc3qOo0wrkvyuhvR2X22F6EXbBRZWTgwAfrXWmUzW7RqqX9R9XWR6KY9ZQO8OAlYZt9CGei15XbXyBRi0aYOBRuWA3vkT8XFoYW2JZ2dDefCHXUGuxhxbor3hVoiQVLZJ6vNTN4wEGWKn1QZaMMHY3D8P7f8f93KaQgg2jXJob9ZTf5Jh1uGvSFPrwh%2BRuEFiu2aDuc5v7vIeI0SW5KhJ7Yimf8mQ8fAM7QNwLLYR3OAbZ%2BSvH0wVJkUn1dm0rw0vRILC7%2BLKLqlzGcSb5ZIrg10hiOaJHg5ns74ps9rHi9giVnUSALt89N2vpZkA75yJf3Q7vJhGZtEJBC8XoI2Y3sWL4zkFwwe%2BrqeEVw%2FZgY1F1BK1AJQkJIWxFTE9hyA%2BEyb70mo710%2BuujyG43e9yrXCHreAgzZ0NOvECMMcvToK%2BbLCufmuWPymxggm%2BGB5wUn2eEu1%2FdZNoKPdlw1hEx0bZ8XwfioijkiwS9d%2B0vGur8BneQSHyVnNj4VO0WrK7WSF4%2FtxxgnFxSKyNfodtXTNVpyn3kwjqe5owY6pAFuaWX3%2FyB0sRcR5ill3p6pKfQHmpF02sPtRarY5SAfdqwr4MPohMh7nR212gIRM1V1sFnMfqO1IlU1VbrzxxFmkHYdes53es54d3RPN4WXk%2BuNWw3yK4N0wdcAlZUeKg4WoaKxtQ8alD0U8GkpYzkyNM4a5Xo%2FXyqr1z991QaykNK%2BFqFAGzvn0ljkXZ%2FrzHRQY57e2BSbosgX55pCHMKVqCt9jQ%3D%3D&Expires=1685569384"')
# system('curl -o ref-data/SEAAD/Micro-PVM.rds "https://corpora-data-prod.s3.amazonaws.com/b74ffc69-d983-4e73-a55f-45aef2d7a182/local.rds?AWSAccessKeyId=ASIATLYQ5N5X3ZHE2I4F&Signature=qpGVDqDkWev7Z4uiue%2BsMK0ac2o%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMr%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIHdf5VzJg1FyGovvU0cHWAM9XKZNXW1K%2BLwojfsrFFNNAiAZrA4p%2BbCu54jc3WY7jKvMVwi0bqD9UIMmawCxs8L4vir0Awjz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAEaDDIzMTQyNjg0NjU3NSIMTgX4QAn2507IRb42KsgDySLxEpgxc6jSL2T6%2FRe7WOt52hmanaGBHFHIhrI7nJFQoF%2FaLhOv3X%2FRq%2FNuA%2BQ2HzkbVl0SmwVLfv%2FXM18ac2THqUEfFodykvWulpFB3KnjMv5%2F17RF2BihKDwT2k%2F4KSBE1TaqLjoi90K2rNPAVJH2JdM85SWyRphgdQnFpiF8gE37uvk7ghB8Z%2BMNMHIVSwK%2F3U4TiXN0zOOqUpm2W4gNHW2G2nSade8l1PNtV8KiEHs%2FfHMv42SXiJCqEuPLdusUUsF%2BIUf%2FiAYdqhKTGicQgMAB0Ql0XuuVFmPF3j1M7BmhZSd0VZc6SVYPiPX%2F2lY190r4JgqKXn25XBbaaGqJ7b9JOKKLV61LBNr2X0FXGefK5M6S2EE56l9xqldsufogHfTXWljUCSlr6ikf5t9BRdbAUYRQ6dGBLgS2kIrraTr9F6jskVFC6006e7xPpsCPT%2BNF6%2B1Lmezr7Cpr4A%2FQAgpfO40ie2vfPhCi8IATKr661bG06SuGPy93Pz%2BzW3boZSAy3myWUJDE5c2b%2FWP0OkWR0M5S9sh0W0FIdYy%2BG8w7FgKhIYkH%2FaI%2FSoDPMdXk4nJD%2BGN6lqsEbJNXAVGIVocFBWKxMKeQuaMGOqYBA%2FgL5T6mUs20wCnqgO2MaqwggXOSGBU9aqDNIoFqEPOmsgnskLBQUB5OpJ9n1r0y4W0Do%2BJgcUMxLoZOQTnAff6igSyX6MpsNbEgef8WfF4siBi2IyMbsRRGARXVUGMZkE6WAFi5r%2FZBGgkhL%2FdOhoXtrm%2F2W5%2Bzs5kFWEdDTri1JdRO5lYKFxZE3eeV%2BSegF5FQfgnwNryqOF3fV0LxUzdWb2q%2B4g%3D%3D&Expires=1685569670"')
# system('curl -o ref-data/SEAAD/OPC.rds "https://corpora-data-prod.s3.amazonaws.com/7d763002-8dba-409a-8d71-df1ba61a3cfa/local.rds?AWSAccessKeyId=ASIATLYQ5N5X7ILTHZSG&Signature=rrUhlY3nT5mM8oPm3bNu3%2BKzl6o%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCmk%2BDKAO6TSaqlAMyp8xuGfdAUo6xyoyQX2FL0QymLcwIhAMdmF%2FPA1SQ4zAJg7nmsHlZ%2F9gjda274if1h0dexdS%2BfKvQDCPX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igz9aRgBioj4Xd7Gd70qyAOU13Cq8CtyTVdeynRyHPtB8bMBkLty%2B5dIVx325rmT5yL5D9NmaB%2BOpvuwMpID%2FmQGm3RjCb7B3wnEwu7x%2Bcu%2BhkBVM%2BPBfdBc6jhdqW8X0nvA%2BmhgJF6qA2GtiI%2BuHfapUyMTKupUrbGd%2Fr8wTNjqZj%2B%2BHLs11Mr2vIlIIh0iK2NrMr7BHHpjew25iVfDjoJmvpmnJCI7mLZuvSkXLFhM5SDCMkVscQ%2BX1BrTBV09ImP9uxQKz%2B5AIAUHZaqNxWGrguieesqWFTwsRBhuUyONegDF9nSipNrHpj0G4OlcGdO6pyM2yvcenAfxONM7pVpSl6kZ9SGFGbZU91ggjfSuu3nrfHm%2B11QEuNCn62MPru4il%2FrWAOjy2lYti6j%2B17I2j19m5fw3%2FJcsTuF8f9ZIN0wEMUMMLpWAYFHALV%2FGDVmu60WH2cAKVITeCOkdp7d7r%2BK4KQw7yg8PgPsIJlfTRoBK5rywiUuV9z8VfmioGVjqGkLAN0DL6xSIopTGvOSnlgFIjwf0mRa9OV5mHOb1DsLP5UElXmlA4%2BpRZKhQL%2BKxupPxn1SHZWS6sGEeeVYF0pls%2FRfIPO3u5wlkqv1Iq9wPg4JNd4Mw%2Bse5owY6pAFuZMwAir42WfJK0XtVwhqPRxdzcsNGdFHypq%2BOYwbbOrNRnqxg3A13kjpPaY3%2B118HmaMA2YbqLb9%2BWBDHjyD7iRr2SjX%2BMLcAa0UxnOP%2FFWAdSIv2kvuuyb4W0kuCXmAZNkIPdOwyO%2FCyU3spBg40jDpnlJTFu%2BsFmz3if3yFIO4KL%2Fed3WGxXjqBcdDMAhs5D9sNwZsNg6hIUD9GqFybPHsC3w%3D%3D&Expires=1685569728"')
# system('curl -o ref-data/SEAAD/Sncg.rds "https://corpora-data-prod.s3.amazonaws.com/771c03b2-3623-4ea5-a97b-23dcc075f13f/local.rds?AWSAccessKeyId=ASIATLYQ5N5X7ILTHZSG&Signature=UQDnPc7qlw7CNMhKHOB7H39dZto%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCmk%2BDKAO6TSaqlAMyp8xuGfdAUo6xyoyQX2FL0QymLcwIhAMdmF%2FPA1SQ4zAJg7nmsHlZ%2F9gjda274if1h0dexdS%2BfKvQDCPX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igz9aRgBioj4Xd7Gd70qyAOU13Cq8CtyTVdeynRyHPtB8bMBkLty%2B5dIVx325rmT5yL5D9NmaB%2BOpvuwMpID%2FmQGm3RjCb7B3wnEwu7x%2Bcu%2BhkBVM%2BPBfdBc6jhdqW8X0nvA%2BmhgJF6qA2GtiI%2BuHfapUyMTKupUrbGd%2Fr8wTNjqZj%2B%2BHLs11Mr2vIlIIh0iK2NrMr7BHHpjew25iVfDjoJmvpmnJCI7mLZuvSkXLFhM5SDCMkVscQ%2BX1BrTBV09ImP9uxQKz%2B5AIAUHZaqNxWGrguieesqWFTwsRBhuUyONegDF9nSipNrHpj0G4OlcGdO6pyM2yvcenAfxONM7pVpSl6kZ9SGFGbZU91ggjfSuu3nrfHm%2B11QEuNCn62MPru4il%2FrWAOjy2lYti6j%2B17I2j19m5fw3%2FJcsTuF8f9ZIN0wEMUMMLpWAYFHALV%2FGDVmu60WH2cAKVITeCOkdp7d7r%2BK4KQw7yg8PgPsIJlfTRoBK5rywiUuV9z8VfmioGVjqGkLAN0DL6xSIopTGvOSnlgFIjwf0mRa9OV5mHOb1DsLP5UElXmlA4%2BpRZKhQL%2BKxupPxn1SHZWS6sGEeeVYF0pls%2FRfIPO3u5wlkqv1Iq9wPg4JNd4Mw%2Bse5owY6pAFuZMwAir42WfJK0XtVwhqPRxdzcsNGdFHypq%2BOYwbbOrNRnqxg3A13kjpPaY3%2B118HmaMA2YbqLb9%2BWBDHjyD7iRr2SjX%2BMLcAa0UxnOP%2FFWAdSIv2kvuuyb4W0kuCXmAZNkIPdOwyO%2FCyU3spBg40jDpnlJTFu%2BsFmz3if3yFIO4KL%2Fed3WGxXjqBcdDMAhs5D9sNwZsNg6hIUD9GqFybPHsC3w%3D%3D&Expires=1685569837"')
# system('curl -o ref-data/SEAAD/Lamp5_Lhx6.rds "https://corpora-data-prod.s3.amazonaws.com/be77a165-71a8-4d0c-93f0-f969c213162f/local.rds?AWSAccessKeyId=ASIATLYQ5N5X7ILTHZSG&Signature=LcYts03SYZvIECEgyd%2BN3TUSwEM%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCmk%2BDKAO6TSaqlAMyp8xuGfdAUo6xyoyQX2FL0QymLcwIhAMdmF%2FPA1SQ4zAJg7nmsHlZ%2F9gjda274if1h0dexdS%2BfKvQDCPX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igz9aRgBioj4Xd7Gd70qyAOU13Cq8CtyTVdeynRyHPtB8bMBkLty%2B5dIVx325rmT5yL5D9NmaB%2BOpvuwMpID%2FmQGm3RjCb7B3wnEwu7x%2Bcu%2BhkBVM%2BPBfdBc6jhdqW8X0nvA%2BmhgJF6qA2GtiI%2BuHfapUyMTKupUrbGd%2Fr8wTNjqZj%2B%2BHLs11Mr2vIlIIh0iK2NrMr7BHHpjew25iVfDjoJmvpmnJCI7mLZuvSkXLFhM5SDCMkVscQ%2BX1BrTBV09ImP9uxQKz%2B5AIAUHZaqNxWGrguieesqWFTwsRBhuUyONegDF9nSipNrHpj0G4OlcGdO6pyM2yvcenAfxONM7pVpSl6kZ9SGFGbZU91ggjfSuu3nrfHm%2B11QEuNCn62MPru4il%2FrWAOjy2lYti6j%2B17I2j19m5fw3%2FJcsTuF8f9ZIN0wEMUMMLpWAYFHALV%2FGDVmu60WH2cAKVITeCOkdp7d7r%2BK4KQw7yg8PgPsIJlfTRoBK5rywiUuV9z8VfmioGVjqGkLAN0DL6xSIopTGvOSnlgFIjwf0mRa9OV5mHOb1DsLP5UElXmlA4%2BpRZKhQL%2BKxupPxn1SHZWS6sGEeeVYF0pls%2FRfIPO3u5wlkqv1Iq9wPg4JNd4Mw%2Bse5owY6pAFuZMwAir42WfJK0XtVwhqPRxdzcsNGdFHypq%2BOYwbbOrNRnqxg3A13kjpPaY3%2B118HmaMA2YbqLb9%2BWBDHjyD7iRr2SjX%2BMLcAa0UxnOP%2FFWAdSIv2kvuuyb4W0kuCXmAZNkIPdOwyO%2FCyU3spBg40jDpnlJTFu%2BsFmz3if3yFIO4KL%2Fed3WGxXjqBcdDMAhs5D9sNwZsNg6hIUD9GqFybPHsC3w%3D%3D&Expires=1685569932"')
#
# system('curl -o ref-data/SEAAD/Chandelier.rds "https://corpora-data-prod.s3.amazonaws.com/7b6c7301-7a1d-4096-ba48-7baf7ed585c1/local.rds?AWSAccessKeyId=ASIATLYQ5N5X7ILTHZSG&Signature=g4P3FDJIc9OGUe0LRHwrI54%2FbCc%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCmk%2BDKAO6TSaqlAMyp8xuGfdAUo6xyoyQX2FL0QymLcwIhAMdmF%2FPA1SQ4zAJg7nmsHlZ%2F9gjda274if1h0dexdS%2BfKvQDCPX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igz9aRgBioj4Xd7Gd70qyAOU13Cq8CtyTVdeynRyHPtB8bMBkLty%2B5dIVx325rmT5yL5D9NmaB%2BOpvuwMpID%2FmQGm3RjCb7B3wnEwu7x%2Bcu%2BhkBVM%2BPBfdBc6jhdqW8X0nvA%2BmhgJF6qA2GtiI%2BuHfapUyMTKupUrbGd%2Fr8wTNjqZj%2B%2BHLs11Mr2vIlIIh0iK2NrMr7BHHpjew25iVfDjoJmvpmnJCI7mLZuvSkXLFhM5SDCMkVscQ%2BX1BrTBV09ImP9uxQKz%2B5AIAUHZaqNxWGrguieesqWFTwsRBhuUyONegDF9nSipNrHpj0G4OlcGdO6pyM2yvcenAfxONM7pVpSl6kZ9SGFGbZU91ggjfSuu3nrfHm%2B11QEuNCn62MPru4il%2FrWAOjy2lYti6j%2B17I2j19m5fw3%2FJcsTuF8f9ZIN0wEMUMMLpWAYFHALV%2FGDVmu60WH2cAKVITeCOkdp7d7r%2BK4KQw7yg8PgPsIJlfTRoBK5rywiUuV9z8VfmioGVjqGkLAN0DL6xSIopTGvOSnlgFIjwf0mRa9OV5mHOb1DsLP5UElXmlA4%2BpRZKhQL%2BKxupPxn1SHZWS6sGEeeVYF0pls%2FRfIPO3u5wlkqv1Iq9wPg4JNd4Mw%2Bse5owY6pAFuZMwAir42WfJK0XtVwhqPRxdzcsNGdFHypq%2BOYwbbOrNRnqxg3A13kjpPaY3%2B118HmaMA2YbqLb9%2BWBDHjyD7iRr2SjX%2BMLcAa0UxnOP%2FFWAdSIv2kvuuyb4W0kuCXmAZNkIPdOwyO%2FCyU3spBg40jDpnlJTFu%2BsFmz3if3yFIO4KL%2Fed3WGxXjqBcdDMAhs5D9sNwZsNg6hIUD9GqFybPHsC3w%3D%3D&Expires=1685570014"')
#
# system('curl -o ref-data/SEAAD/Pax6.rds "https://corpora-data-prod.s3.amazonaws.com/0c444642-5403-429b-aa2d-f4f9fed8f25d/local.rds?AWSAccessKeyId=ASIATLYQ5N5X572JEX7C&Signature=IEqFe7BK08UzOkfz7k2i5ZtECw8%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMr%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQCSLYD6c9GJDJ462R3fyKA%2BtqLAfd3ZLYdTs2l8Sm8OqAIhAOufFKbvnALqcS51UPJRams4iUvBfeoo4DEazM%2BCELe8KvQDCPP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1IgzVqIuk562TH8S9imEqyAN32f3C0LsQrgnOZHSbd%2F35aLvFO88avQ3fpI%2B9aFN4d4cQhj4%2FAXOaGvo1pTcjmSZvzc3qOo0wrkvyuhvR2X22F6EXbBRZWTgwAfrXWmUzW7RqqX9R9XWR6KY9ZQO8OAlYZt9CGei15XbXyBRi0aYOBRuWA3vkT8XFoYW2JZ2dDefCHXUGuxhxbor3hVoiQVLZJ6vNTN4wEGWKn1QZaMMHY3D8P7f8f93KaQgg2jXJob9ZTf5Jh1uGvSFPrwh%2BRuEFiu2aDuc5v7vIeI0SW5KhJ7Yimf8mQ8fAM7QNwLLYR3OAbZ%2BSvH0wVJkUn1dm0rw0vRILC7%2BLKLqlzGcSb5ZIrg10hiOaJHg5ns74ps9rHi9giVnUSALt89N2vpZkA75yJf3Q7vJhGZtEJBC8XoI2Y3sWL4zkFwwe%2BrqeEVw%2FZgY1F1BK1AJQkJIWxFTE9hyA%2BEyb70mo710%2BuujyG43e9yrXCHreAgzZ0NOvECMMcvToK%2BbLCufmuWPymxggm%2BGB5wUn2eEu1%2FdZNoKPdlw1hEx0bZ8XwfioijkiwS9d%2B0vGur8BneQSHyVnNj4VO0WrK7WSF4%2FtxxgnFxSKyNfodtXTNVpyn3kwjqe5owY6pAFuaWX3%2FyB0sRcR5ill3p6pKfQHmpF02sPtRarY5SAfdqwr4MPohMh7nR212gIRM1V1sFnMfqO1IlU1VbrzxxFmkHYdes53es54d3RPN4WXk%2BuNWw3yK4N0wdcAlZUeKg4WoaKxtQ8alD0U8GkpYzkyNM4a5Xo%2FXyqr1z991QaykNK%2BFqFAGzvn0ljkXZ%2FrzHRQY57e2BSbosgX55pCHMKVqCt9jQ%3D%3D&Expires=1685570456"')
# system('curl -o ref-data/SEAAD/VLMC.rds "https://corpora-data-prod.s3.amazonaws.com/97b53668-3e54-4436-a1ca-b18e93183323/local.rds?AWSAccessKeyId=ASIATLYQ5N5X3ZHE2I4F&Signature=XFqSkm96VX09CqQFN9yaTKGUcr0%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMr%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIHdf5VzJg1FyGovvU0cHWAM9XKZNXW1K%2BLwojfsrFFNNAiAZrA4p%2BbCu54jc3WY7jKvMVwi0bqD9UIMmawCxs8L4vir0Awjz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAEaDDIzMTQyNjg0NjU3NSIMTgX4QAn2507IRb42KsgDySLxEpgxc6jSL2T6%2FRe7WOt52hmanaGBHFHIhrI7nJFQoF%2FaLhOv3X%2FRq%2FNuA%2BQ2HzkbVl0SmwVLfv%2FXM18ac2THqUEfFodykvWulpFB3KnjMv5%2F17RF2BihKDwT2k%2F4KSBE1TaqLjoi90K2rNPAVJH2JdM85SWyRphgdQnFpiF8gE37uvk7ghB8Z%2BMNMHIVSwK%2F3U4TiXN0zOOqUpm2W4gNHW2G2nSade8l1PNtV8KiEHs%2FfHMv42SXiJCqEuPLdusUUsF%2BIUf%2FiAYdqhKTGicQgMAB0Ql0XuuVFmPF3j1M7BmhZSd0VZc6SVYPiPX%2F2lY190r4JgqKXn25XBbaaGqJ7b9JOKKLV61LBNr2X0FXGefK5M6S2EE56l9xqldsufogHfTXWljUCSlr6ikf5t9BRdbAUYRQ6dGBLgS2kIrraTr9F6jskVFC6006e7xPpsCPT%2BNF6%2B1Lmezr7Cpr4A%2FQAgpfO40ie2vfPhCi8IATKr661bG06SuGPy93Pz%2BzW3boZSAy3myWUJDE5c2b%2FWP0OkWR0M5S9sh0W0FIdYy%2BG8w7FgKhIYkH%2FaI%2FSoDPMdXk4nJD%2BGN6lqsEbJNXAVGIVocFBWKxMKeQuaMGOqYBA%2FgL5T6mUs20wCnqgO2MaqwggXOSGBU9aqDNIoFqEPOmsgnskLBQUB5OpJ9n1r0y4W0Do%2BJgcUMxLoZOQTnAff6igSyX6MpsNbEgef8WfF4siBi2IyMbsRRGARXVUGMZkE6WAFi5r%2FZBGgkhL%2FdOhoXtrm%2F2W5%2Bzs5kFWEdDTri1JdRO5lYKFxZE3eeV%2BSegF5FQfgnwNryqOF3fV0LxUzdWb2q%2B4g%3D%3D&Expires=1685570525"')
# system('curl -o ref-data/SEAAD/endo.rds "https://corpora-data-prod.s3.amazonaws.com/b9968dc3-5afb-4f26-9b60-524096ede646/local.rds?AWSAccessKeyId=ASIATLYQ5N5XTJULZ6FX&Signature=cDiJR%2Fv42T%2B167ecOPIQ2XRF830%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIHANI3tq2E2kkqgt9rN9bFItFqqLjG%2F%2BukVNlqm21qr0AiA6wTlG6NhNqvyVqxQ9AEETFBhyKfNvzxei0yWNJ7ah1Sr0Awj0%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAEaDDIzMTQyNjg0NjU3NSIMJlE6hpbBIxFtFYkyKsgDtquTcniXhHXdRKg%2FD6oVy15CNAm9ad16qfpdPrmuv%2B8cnhVZGA0%2FVeWvza49oG4Wgj%2B%2BBt5T7wsKWdJbPbrfIDIkb2C5c93B0LsRAxDemoGOtJ7uyD3Gz%2FIqyD7pspOxGtuQdIaeMMtaTOXsumn3Q9wR%2FrmofTnbbx3Re3oeKmemb27QYTGHJxFKflx%2B%2BJuQLn6OL%2Bz1GeNRtaCq9g2iLfkV7lWzqHBM5MiqpuMRIF1bPSbjDFVAkK8szfmHVC9JO%2BtBJsOiAjauXTYIGLYWBBYl%2B%2FOhjrvwIxOkyzQDcFPXjO2Xi0QBdSmSWjWswTz8dZEKggap2vOqea4mnRPNU3ZLAHhWPoZfokcvY4MImqpVWza%2BJO5BRsutgUCiwgMWqE%2FM5PhIigACj1cZIo%2F1pRvhTK7wcgCxotRpZ1EqPn4ctgLH4tNUbpwulI4S%2Fbpmr%2FpHImpJ9hLYgQv5XiSiKAumGelOnqNdW3qVFk5XxR%2Ft6aSTCxQs2zmbCyZOk9Kd6lOai4LgaR1tHUmrH%2FB1DlrMh%2B7g3Bpe7bRHHippTRx7c%2BXqSxF7317XOhXkyeHWIObnOtmFBoiqlybKNq72RwuGwIcsdBCsMMbGuaMGOqYBo0dUCzf6veK4K7I%2FZPdhwR1nglVfOnGxCOwQd3na4ks%2F4gFRPd%2FLoEndEgSaip%2FSA%2FiNfIngGL34yiadYo9eJ9xXxxK%2FUJ0eP%2BCAdn44PTCw893lg5fXWhVU5BoK0cwbzybnZKcctlrL4XxVVcrgiy0Kv8NO5y6Ck9Qy4S0LRLvhUm0k45qSnwqQdOpMQsv0u7yHmXxhwBn9AfIOcom8eTO4I240Lg%3D%3D&Expires=1685570583"')
# system('curl -o ref-data/SEAAD/Sst-Chodl.rds "https://corpora-data-prod.s3.amazonaws.com/7365e4e1-4711-494b-8171-fda7f5611872/local.rds?AWSAccessKeyId=ASIATLYQ5N5XV2SDXEL5&Signature=um3v%2BiaeVeJTFZrzxPXgUwb2RQc%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEMv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCID3h3jAtF8yhD5RZrNPzmbhbWTfErrDZfhoSWWkeIwCcAiEAzAMU0FhETPpyBAebZTV%2BjbCcGav7T7%2BeYO%2B38iZIbMkq9AMI8%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgwyMzE0MjY4NDY1NzUiDPATOnyPkr085qmBzCrIA3Bwhue0%2FaxyfGFqOE1Y2IFqZiOYQ65p73TFMHsr%2FdnCQwehzEA%2BsWxCpieiyWCZBN61%2FjxE%2FT3xhVGIj16QyGQqiLAkuWPkvNmp9o60urX5po59WgU1CyX9Gp8Kx6rra3rhE5ckzmUlFU4vkiClPZEc6baQKVxqPIuLYynFlnaDXwtWZdvz35nGVkGVLSG53yyTFhb79GSFaMuTDosZS%2BZ686dtRUaxGCTbRm%2FNBhIU1LEmp2WkB9qpo5Tj3lRC7OEXSN0%2FwGEWhCMgrhNdLFq9Xi342ZqjO7SkcSSQE85u5EmY2kgl887tTzmV4QCQx4Bm6AcSZBes%2FneWWSBboaZRHgs%2Fnp29n%2F8BCzp%2BMpeE8pgwHNwV%2BzS%2Bm5d3u72%2FOBqNjtEveuLU5cGhLnxBjs5rNOAK%2Bu0XwC%2BsJbz3Ah42vsVqPjW%2FoXc6fUBZsOF6C15BCXmIqSKeTQHLXCdNTHgC0nJfumMUJio99GBhn9uIGB5x04hcX9gp6YRoUk2IEI2Oj3XaGC7nWZLNtPxUGSBJ5ne8suEtpFlaRZMIRUWbhns2uf7T4rbLVWIb1j%2BJdbsww6FgSUmLJ8aB0WPTqSR2SSsjj3BTkzCeqrmjBjqlAevXbGwDpQ03KR6r02VfrEuPhw%2B98UqZuKFh9%2BX13Ootdr4XGWMkK0HRkqMlij2zsla%2BkQL4TZjgANz3n6Q4nAAry7dKfDtZsm06Luc6U1sEJWcbMG19Q08BTocCFr5uEKkXqKPYklCoApXjldV3Sm6odP8A4EvAsUbB9%2BP6bOlpkDlPwCBrVIiqDdVRgVD9bCGPErJ77GYAy2NJy84lVuOyw2q%2B7w%3D%3D&Expires=1685570644"')
#

#
#produce a 50k proport. cells datasets
# cells<-readRDS('ref-data/SEAAD/endo.rds')
# head(cells@meta.data)
# rm(cells)

brains100k<-Reduce(function(x,y)merge(x,y,merge.dr=c('scVI','umap')),lapply(file.path('ref-data/SEAAD',list.files('ref-data/SEAAD/',pattern = '.rds')),function(fp){
  obj<-readRDS(fp)
  #subset div by 10
  mtd<-data.table(obj@meta.data,keep.rownames = 'bc')
  mtd[,to_keep:=bc%in%sample(bc,round(.N/10))]
  objf<-obj[,mtd[(to_keep)]$bc]
  return(objf)
}))
message('proportional 100k cells object creation finished')

table(brains100k$cell_type)


#save
saveRDS(brains100k,fp(out,'SEAAD_brains_MTG_100k.rds'))


#rm(brains100k)

#20k subset : take at least500 by cell_type + 20%of , max 2k cell type

brains20k<-Reduce(function(x,y)merge(x,y,merge.dr=c('scVI','umap')),lapply(file.path('ref-data/SEAAD',list.files('ref-data/SEAAD/',pattern = '.rds')),function(fp){
  obj<-readRDS(fp)
  #subset div by 10
  mtd<-data.table(obj@meta.data,keep.rownames = 'bc')
  mtd[,to_keep:=bc%in%head(sample(bc,ifelse(.N*0.2>500,round(.N*0.2),500)),2000)]
  
  objf<-obj[,mtd[(to_keep)]$bc]
  return(objf)
}))

#DimPlot(brains20k,group.by = 'cell_type')
message('500cells 2kcells max 20% downsampling object creation finished')

#save
saveRDS(brains20k,fp(out,'SEAAD_brains_MTG_20k.rds'))

message('finish!')


#see the results
brain100k<-readRDS('outputs/01-SEAAD_data/SEAAD_brains_MTG_100k.rds')

DimPlot(brain100k,group.by='cell_type',label=T,reduction = 'umap') #not the same umap, need to redo clustering
head(brain100k@meta.data)
mtd<-data.table(brain100k@meta.data,keep.rownames = 'bc')
table(unique(mtd,by='donor_id')$disease)
# dementia   normal 
#       42       47 


brain20k<-readRDS(fp(out,'SEAAD_brains_MTG_20k.rds'))

brain20k<-NormalizeData(brain20k)

brain20k<-FindVariableFeatures(brain20k)

brain20k<-ScaleData(brain20k)

brain20k<-RunPCA(brain20k)

brain20k<-RunUMAP(brain20k,dims = 1:30)
brain20k<-FindNeighbors(brain20k,dims = 1:30)
brain20k<-FindClusters(brain20k,resolution = 0.2)

wrap_plots(lapply(DimPlot(brain20k,group.by=c('cell_type','seurat_clusters'),label = T,
                          label.size = 3,combine = F),function(p)p+NoLegend()))



#only on 3' data
brain20kf<-subset(brain20k,assay!='10x multiome')
brain20kf


#stats
table(brain20kf$donor_id)
mtd<-data.table(brain20kf@meta.data,keep.rownames = 'cell_id')
table(unique(mtd,by='donor_id')$disease)
# dementia   normal 
#       42       47 

#save object
saveRDS(brain20kf,fp(out,'SEAAD_brains_MTG_34k.rds'))


#save matrix and metadata for exercise

library(Matrix)
out1<-fp(out,'Young_coding_exercise')
dir.create(out1)
# save sparse matrix
out2<-fp(out1,'SEA_AD34k_filtered_feature_bc_matrix')
dir.create(out2)

writeMM(obj = brain20kf@assays$RNA@counts, file=fp(out2,"matrix.mtx"))
#system: gzip matrix.mtx 

# save genes and cells names
fwrite(list(rownames(brain20kf@assays$RNA@counts)), file = fp(out2,"features.tsv.gz"))
fwrite(x = list(colnames(brain20kf@assays$RNA@counts)), file = fp(out2,"barcodes.tsv.gz"))


#metadata without cell type info
head(brain20kf[[]])
fwrite(data.table(brain20kf@meta.data,keep.rownames = 'cell_id')[,.(cell_id,donor_id,tissue,disease,sex,self_reported_ethnicity,`Age at death`,`Years of education`,`APOE4 status`)],
       fp(out1,'SEA_AD34k_metadata.csv.gz'))

#test
mat<-Read10X(out2,gene.column = 1)
dim(mat)
head(mat[,1:10])
colnames(mat)
mtd<-fread(fp(out1,'SEA_AD34k_metadata.csv.gz'))
mtd     
#OK
